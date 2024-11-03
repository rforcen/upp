#include <Core/Core.h>
#include <CtrlLib/CtrlLib.h>
#include <GridCtrl/GridCtrl.h>

#include <algorithm>

#include "../numCpp/numCpp.h"
#include "../numCpp/Timer.h"

using namespace Upp;
using namespace std;

using Real = double;

class NCWin : public TopWindow {
 private:
  class NCGrid : public GridCtrl {
   private:
    NC<Real>* nc = nullptr;
    VInt slice;

   public:
    NCGrid() { Sorting(false); }

    void setNC(NC<Real>& nc, VInt slice) {
      auto rc = GetRect();
      SetRect(0, 0, 0, 0);  // set rect to 0 to force refresh

      this->nc = &nc;
      this->slice = slice;

      Disp();

      SetRect(rc);
    }

    void Disp() {
      Clear(true);

      if (nc != nullptr) {
        for (int c = 0; c <= nc->dim(0); c++)
          AddColumn(c == 0 ? "#" : AsString(c)).Width(120);

        auto slc = (slice.size() > 0) ? nc->slice(slice) : *nc;

        if (slc.ndims == 2) {
          for (int r = 0; r < slc.dim(1); r++) {
            for (int c = 0; c <= slc.dim(0); c++)
              Set(r, c,
                  c == 0 ? AsString(r + 1) : Format("%.3g", slc.at(r, c - 1)));
          }
        } else {
          for (int c = 0; c <= slc.dim0; c++)
            Set(0, c, c == 0 ? AsString(1) : Format("%.3g", slc.at(c - 1)));
        }

        SetColsMin(100);
      }
    }
  };

  // controls
  StatusBar sb;
  ToolBar tb;
  Label nItems;                     // nitems on 'a' matrix
  Vector<EditIntSpin*> dims;        // dims spins
  NCGrid ncGrid;                    // data grid
  EditString sDims;                 // csv dims editor
  EditDouble eConst;                // const to operate with 'a' matrix
  FrameLeft<ParentCtrl> leftframe;  // buttons frame
  ParentCtrl container;             // & its container
  Splitter splitter;
  EditIntSpin nDims;

  Vector<Button*> vbutt;
  // buttStr & enum must be on sync
  Vector<String> buttStr = {"b=a",  "a=b",  "rand", "a+b", "a-b",   "a*b",
                            "a/b",  "det",  "inv",  "dot", "saveA", "saveB",
                            "zero", "ones", "fill", "sum", "norm",  "minMax",
                            "+",    "-",    "*",    "/"};
  typedef enum {
    A_BEQA,
    A_AEQB,
    A_RAND,
    A_APLUSB,
    A_AMINUSB,
    A_APRODB,
    A_ADIVB,
    A_DET,
    A_INV,
    A_DOT,
    A_SAVEA,
    A_SAVEB,
    A_ZERO,
    A_ONES,
    A_FILL,
    A_SUM,
    A_NORM,
    A_MINMAX,
    A_PLUS,
    A_MINUS,
    A_PROD,
    A_DIV
  } Actions;

  NC<Real> a, b;  // working matrix
  VInt slice;     // current slice to display

 public:
  NCWin() {
    Title("NumCpp").Sizeable().CenterScreen().SetRect(Size(1600, 1200));
    leftFrameCook();

    // a = NC<Real>::random(10, 20, 3, 10, 10);  // NC init
    generateRandom(5);

    splitter << leftframe << ncGrid;  // split leftframe / grid
    splitter.SetPos(1300);
    Add(splitter.SizePos());

    sb = "NC";
    AddFrame(sb);  // status bar
    AddFrame(tb);  // tool bar

    refreshGridTB();
  }

  ~NCWin() {
    clearDims();
    for (auto& v : vbutt) delete v;
  }

 private:
  const void refreshGridTB() {  // call everytime 'a' changes
    slice.resize((a.ndims > 2) ? a.ndims - 2 : 0, 0);
    ncGrid.setNC(a, slice);  // grid disp

    toolBarPopulate();
  }
  const void openNPY() {  // from numpy file
    {
      FileSel fs;
      if (fs.Type("numpy files", "*.npy").ExecuteOpen("open npy file")) {
        a.load(fs.Get().ToStd());

        refreshGridTB();
      }
    }
  }
  const void saveNPY() {  // to numpy file
    {
      FileSel fs;
      if (fs.Type("numpy files", "*.npy").ExecuteSaveAs("save npy file"))
        a.save(fs.Get().ToStd());
    }
  }
  void dispGrid() { ncGrid.setNC(a, slice); }

  Real getConst() { return (Real)eConst.GetData(); }
  void doAction(String action) {
    bool upSB = true;
    Timer ts;

    splitter.SetRect(0, 0, 0, 0);

    switch (FindIndex(buttStr, action)) {
      case A_BEQA:
        b = a;
        break;
      case A_AEQB:
        a = b;
        dispGrid();
        break;
      case A_RAND:
        a.rand();
        dispGrid();
        break;
      case A_APLUSB:
        a += b;
        dispGrid();
        break;
      case A_AMINUSB:
        a -= b;
        dispGrid();
        break;
      case A_APRODB:
        a *= b;
        dispGrid();
        break;
      case A_ADIVB:
        a /= b;
        dispGrid();
        break;
      case A_DET:
        a = a.det();
        refreshGridTB();
        break;
      case A_INV:
        a = a.inv();
        dispGrid();
        break;
      case A_DOT:
        a = a.dot(b);
        refreshGridTB();
        break;
      case A_SAVEA:
        a.save("a.npy");
        break;
      case A_SAVEB:
        b.save("b.npy");
        break;
      case A_ZERO:
        a.zeros();
        dispGrid();
        break;
      case A_ONES:
        a.fill(1);
        dispGrid();
        break;
      case A_FILL:
        a.fill(getConst());
        dispGrid();
        break;
      case A_SUM:
        sb = Format("Sum:%g", a.sum());
        upSB = false;
        break;
      case A_NORM:
        a.norm();
        dispGrid();
        break;
      case A_MINMAX: {
        auto mm = a.minMax();
        sb = Format("min max:<%g,%g>", mm.first, mm.second);
        upSB = false;
      } break;
      case A_PLUS:
        a += getConst();
        dispGrid();
        break;
      case A_MINUS:
        a -= getConst();
        dispGrid();
        break;
      case A_PROD:
        a *= getConst();
        dispGrid();
        break;
      case A_DIV:
        a /= getConst();
        dispGrid();
        break;
    };

    splitter.SizePos();

    if (upSB) sb = Format("lap for %s: %d ms", action, (int)ts.lap());
  }
  void leftFrameCook() {  // in leftframe

    int i = 0, hoff = 10, bw = 100, bh = 32, ncol = 2, c, r;

    for (auto& bs : buttStr) {  // button grid
      Button* btn = new Button;

      btn->SetLabel(bs);

      c = bw * (i % ncol);
      r = bh * (i / ncol) + hoff;
      btn->SetRect(c, r, bw, bh);

      *btn << [this, btn] { doAction(btn->GetLabel()); };

      container << *btn;
      vbutt << btn;

      i++;
    }

    // constant to operate w/'a' matrix
    i++;
    eConst.SetRect(0, bh * (i / ncol) + hoff, 200, 32);
    eConst.SetData(1);
    eConst.Tip("constant used in +-*/, fill");
    container << eConst;

    leftframe.Add(container.SizePos());
  }
  void redispGrid(Ctrl& c) {  // force redisplay
    c.Hide();
    c.Show();
  }
  void generateRandom(int nd) {
    VInt dm;
    for (int i = 0; i < nd; i++) dm.push_back(1 + random() % 10);
    if (dm.size() > 1) dm.back() = dm[dm.size() - 2];  // force qudratic

    a = NC<Real>::random(dm);

    refreshGridTB();
  }
  void toolBarPopulate() {  // populate controls in ToolBar tb
    tb.Clear();

    auto ad = pickFirst_2(a);  // first ndims-2 dimensions

    tb.Add("open", CtrlImg::open(), [this] {
        openNPY();
      }).Key(K_CTRL_O);  // open btn
    tb.Add("save", CtrlImg::save(), [this] {
        saveNPY();
      }).Key(K_CTRL_S);  // save btn

    nDims.MinMax(1, 12);  // n dimensions
    nDims.SetRectX(0, 90);
    nDims.SetData((int)a.ndims);
    nDims.Tip("number of dimensions");
    nDims.WhenAction = [this] {  // generate random matrix with nDims dimensions
      generateRandom((int)nDims.GetData());
    };
    tb.Add(nDims, 110);

    sDims.Tip("Dimensions");  // dims editor
    sDims.SetRectX(0, 270);
    sDims.SetText(v2s(a.dims));  // set current dimensions
    sDims.WhenEnter << [this] {
      VInt dm = String2VInt(sDims.GetData());
      if (dm.size() != 0) {
        a = NC<Real>::random(dm);
        refreshGridTB();
      }
    };
    tb.Add(sDims, 270);

    nItems.SetLabel("# items:" + AsString(a.size));  // label
    tb.Add(nItems, 270);

    clearDims();
    int nd = 0;
    for (auto& d : ad) {  // array of spins with selected slice
      EditIntSpin* sp = new EditIntSpin;

      sp->MinMax(0, d - 1);
      sp->SetData(0);
      sp->Tip(Format("dim%d(0..%d)", nd++, (int)(d - 1)));

      sp->WhenAction = [this] {  // Int Spin action
        slice.clear();
        for (auto& d : dims) slice.push_back(*d);

        ncGrid.setNC(a, slice);

        sb = Format("[%s]", v2s(slice));
      };

      sp->SetRectX(0, 90);
      tb.Add(*sp, 90);

      dims.Add(sp);
    }

    redispGrid(tb);
  }
  void clearDims() {
    for (auto& d : dims) delete d;
    dims.Clear();
  }
  String v2s(VInt v) {
    String s;
    for (auto& si : v) s += AsString(si) + ",";
    if (s.EndsWith(",")) s.Trim(s.GetLength() - 1);
    return s;
  }
  VInt String2VInt(const String& str) {
    VInt v;

    regex re(R"(\d+)");
    string s = str.ToStd();

    for (auto it = sregex_token_iterator(s.begin(), s.end(), re);
         it != sregex_token_iterator(); it++)
      v.push_back(stol(*it));

    return v;
  }
  VInt pickFirst_2(NC<Real>& a) {  // pick first n-2 items
    VInt ad;                       // first n-2 items
    if (a.ndims > 2) ad.insert(ad.begin(), a.dims.begin(), a.dims.end() - 2);
    return ad;
  }
};

GUI_APP_MAIN {
  srand(time(nullptr));
  NCWin().Run();
}

