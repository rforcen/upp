/*
    Mandelbrot Fractals explorer
*/
#include <boost/multiprecision/cpp_bin_float.hpp>

#include <CtrlLib/CtrlLib.h>
#include <GridCtrl/GridCtrl.h>

#include "mandel/mandel.h"

using namespace Upp;
using namespace std;

// fp formats
using f32 = float;
using f64 = double;
using f128 = long double;
using f256 =
    boost::multiprecision::number<boost::multiprecision::cpp_bin_float<256>>;

typedef f256 Real;
typedef complex<Real> Cmplx;

class FracWin : public TopWindow {
  typedef enum { pf32, pf64, pf128, pf256 } PrecType;

  struct FractalGeo {  // center, range of fractal in max precission available (f256)
    PrecType pt = pf64;
    complex<f256> center, range;
    int iters = 600;

    FractalGeo() { init(); }

    void init() {
      pt = pf64;
      iters = 600;
      center = complex<f256>(0.5, 0);
      range = complex<f256>(-2, 2);
    }
  } fg;

  struct MySplitter : public Splitter  // splitter that refreshes main win
  {
    FracWin *frWin;

    MySplitter(FracWin *frWin) : frWin(frWin) {}

    virtual void MouseMove(Point p, dword keyflags) override {
      Splitter::MouseMove(p, keyflags);

      frWin->repaint();
    }
  };

  struct MyCtrlImage : public ImageCtrl { // fractal image
    FracWin *frWin;

    MyCtrlImage(FracWin *frWin) : frWin(frWin) {}

    virtual void RightDown(Point p, dword keyflags) override {  // zoom out
      if (keyflags & K_CTRL) {
        frWin->fg.range /= 0.8;
        frWin->changed = true;
        frWin->repaint();
      } else {
        MenuBar::Execute(
            this,
            [=](Bar &bar) {
              const vector<pair<int, String>> resols = {
                  {K_ALT_N, "new"},
                  {K_F2, "save session"},
                  {K_F3, "open saved session"},
                  {K_F4, "toggle mosaic"},
                  {K_F5, "save image"},
                  {K_F6, "save all mosaic images"},
                  {K_ALT_ENTER, "save to selected"}};

              for (auto &r : resols)
                bar.Add(GetKeyDesc(r.first) + " " + r.second,
                        [=] { frWin->Key(r.first, 1); });
            },
            p + frWin->GetRect().TopLeft());
      }
    }

    virtual void LeftDown(Point p, dword keyflags) override {  // zoom to point
      ImageCtrl::LeftDown(p, keyflags);

      frWin->reCalculateCenterRange(p);
    }
  };

  struct MyGridCtrl : public GridCtrl { // left grid with fractal visited points
    FracWin *frWin;

    MyGridCtrl(FracWin *frWin) : frWin(frWin) {}

    virtual bool Key(dword key, int n) {
      switch (key) {
        case K_DELETE:
          if (GetCurrentRow() >= 0) {
            int cr = GetCurrentRow();
            Remove(cr);
          }
          return true;
      }
      return frWin->Key(key, n);  // GridCtrl::Key(key, n);
    }
  };

  struct MosaicCtrl : public Ctrl { // mosaic view of all visited points
    FracWin *frWin;
    ScrollBar sb;

    MosaicCtrl(FracWin *frWin) : frWin(frWin) {
      AddFrame(sb);
      sb.WhenScroll = [=] { Refresh(); };
    }

    void redispSB() {
      int h = GetSize().cy, w = GetSize().cx;
      int grSize = frWin->grSize;

      sb.SetLine(grSize);

      int nc = frWin->visGrid.GetCount() / (w / grSize);

      sb.SetTotal(grSize * nc);
      sb.SetPage(h);
    }

    void LeftDown(Point p, dword flags) override {
      int grSize = frWin->grSize;
      int nImg = p.x / grSize + (p.y / grSize) * (GetSize().cx / grSize) +
                 sb / grSize;  // offset in visGrid;

      if (nImg < frWin->visGrid.GetCount()) {
        frWin->fg = ValueTo<FractalGeo>(frWin->visGrid.Get(nImg, frWin->colFG));
        frWin->toggleDisplayMode();  // toogle to Single
      }
      Ctrl::LeftDown(p, flags);
    }

    void Paint(Draw &dw) override {
      redispSB();

      int grSize = frWin->grSize;
      int w = GetSize().cx, h = GetSize().cy, nc = w / grSize;

      int voff = sb / grSize;  // offset in vImage

      dw.DrawRect(GetRect(), White);
      for (int r = 0; r + voff < frWin->visGrid.GetCount();
           r++) {  // redisp vImage
        int x = r % nc, y = r / nc;

        Image tn =
            ValueTo<Image>(frWin->visGrid.Get(r + voff, frWin->colImage));

        dw.DrawImage(x * grSize, y * grSize, tn);
      }
    }
  };

 private:
  StatusBar sb;
  MySplitter spl;
  MyCtrlImage fracDisp;
  MyGridCtrl visGrid;

  MosaicCtrl mosaic;

  EditInt edIters;  // sb controls
  Label lsb, l0, l1, l2;
  EditIntSpin res;
  DropList precs;

  TimeStop ts;  // lap
  double lap;

  volatile bool resizing = false;
  bool changed = true;

  Image imgFrac;
  typedef enum { Single, Mosaic } DispMode;
  DispMode dmode = Single;

  // grid
  const int grSize = 256, colIC = 0, colFG = 1, colImage = 2;

  typedef FracWin CLASSNAME;

 public:
  FracWin() : spl(this), fracDisp(this), visGrid(this), mosaic(this) {
    Title("Mandelbrot Fractals")
        .Sizeable()
        .Zoomable()
        .CenterScreen()
        .BackPaint()
        .SetRect(Size(6 * grSize, 5 * grSize));

    sb = " ";
    AddFrame(sb);  // status bar w/ iter edit
    int hsb = sb.GetSize().cy, xoff = 0;

    auto addLbl = [&](Label &l, String s) {
      s = " " + s + " ";
      l = s;
      l.SetInk(Red());
      int h = GetTextSize(s, StdFont()).cx;

      l.SetRect(xoff, 0, h, hsb);
      xoff += h;
      sb << l;
    };

    addLbl(l0, "iters");

    edIters <<= fg.iters;  // iters
    edIters.SetRect(xoff, 0, 100, hsb);
    xoff += 100;
    sb << edIters;

    addLbl(l2, "prec");

    // precs
    for (auto i : {"f32", "f64", "f128", "f256"}) precs.Add(i);
    precs.Tip("floating point precission");
    precs.SetRect(xoff, 0, 110, hsb);
    xoff += 110;
    precs.SetIndex(1);  // f64
    precs << [=] {
      fg.pt = (PrecType)precs.GetIndex();
      repaint();
    };
    sb << precs;

    addLbl(l1, "resolution");

    res.MinMax(0, 32);  // res: 1..32k
    res.SetInc(1);      // Set increment value
    res.SetRect(xoff, 0, 90, hsb);
    xoff += 90;
    res.Tip("fractal image saved resolution in K, 0 current image size");
    res <<= 0;  // default = current window size
    sb << res;

    addLbl(lsb, String("_", 200));

    // main win controls
    mosaic.Hide();
    Add(mosaic.SizePos());

    spl << visGrid << fracDisp;  // splitter ( visGrid, fracDisp )
    spl.Horz().SetPos(2000);
    Add(spl.SizePos());

    visGrid.AddColumn("fractal");  // grid, show only image
    visGrid.AddColumn("FG-CenterRange");
    visGrid.AddColumn("Image");
    visGrid.HideColumn(2);
    visGrid.HideColumn(3);

    visGrid.WhenCursor = [=] {  // get stored  Cmplx as Value
      if (visGrid.GetCursor() >= 0) {
        fg = ValueTo<FractalGeo>(visGrid.Get(visGrid.GetCursor(), colFG));

        repaint();
      }
    };
  }

 private:
  Image genMandel(int w, int h, int iters, Cmplx center,
                  Cmplx range)  // Mandelbrot interface
  {
    vector<u32> vimg;

    ts.Reset();
    {
      switch (fg.pt) {
        case pf32: {
          vimg = Mandelbrot<f32>::genMandelbrot(
              w, h, iters, complex<f32>((f32)center.real(), (f32)center.imag()),
              complex<f32>((f32)range.real(), (f32)range.imag()));
        } break;
        case pf64:
          vimg = Mandelbrot<f64>::genMandelbrot(
              w, h, iters, complex<f64>((f64)center.real(), (f64)center.imag()),
              complex<f64>((f64)range.real(), (f64)range.imag()));
          break;
        case pf128:
          vimg = Mandelbrot<f128>::genMandelbrot(
              w, h, iters,
              complex<f128>((f128)center.real(), (f128)center.imag()),
              complex<f128>((f128)range.real(), (f128)range.imag()));
          break;
        case pf256:
          vimg = Mandelbrot<f256>::genMandelbrot(w, h, iters, center, range);
          break;
      }
    }
    lap = ts.Elapsed() / 1000;

    ImageBuffer ib(w, h);
    memcpy(ib, vimg.data(), vimg.size() * sizeof(u32));

    return (Image)ib;
  }

  virtual void Layout() override  // when resized
  {
    TopWindow::Layout();

    if (!resizing) {
      resizing = true;
      PostCallback(THISBACK(repaint));
    }
  }

  void toggleDisplayMode() {  // Single / Mosaic
    if (dmode == Single) {
      mosaic.Show();
      mosaic.Refresh();

      fracDisp.Hide();
      visGrid.Hide();
      spl.Hide();

      dmode = Mosaic;
    } else {
      mosaic.Hide();

      fracDisp.Show();
      visGrid.Show();
      spl.Show();

      dmode = Single;
    }
    repaint();
  }

  virtual bool Key(dword key, int count) override  // custom key event
  {
    Real deltaMove = abs(fg.range) / 30.0;
    int w, h;
    auto setGeo = [&] {
      if (~res == 0)  // current window size
      {
        w = fracDisp.GetSize().cx;
        h = fracDisp.GetSize().cy;
      } else  // res value * K
      {
        w = h = (int)~res * 1024;
      }
    };

    switch (key) {
      case K_ESCAPE:  // quit
        Close();
        return true;

      case K_F2: {  // save session
        FileSel fs;
        if (fs.Type("save session", "*.fs").ExecuteSaveAs()) {
          FileOut f(fs.Get());
          if (f)
            for (int r = 0; r < visGrid.GetCount(); r++) {
              auto _fg = ValueTo<FractalGeo>(visGrid.Get(r, colFG));
              f.Put(&_fg, sizeof(_fg));
            }
        }
      }
        return true;

      case K_F3: {  // open session
        FileSel fs;
        if (fs.Type("open session", "*.fs").ExecuteOpen()) {
          FileIn f(fs.Get());
          if (f) {
            visGrid.Clear();  // clear

            while (!f.IsEof()) {
              f.Get(&fg, sizeof(fg));
              addImgGrid();
            }
          }
        }
      }
        return true;

      case K_F4:  // toggle display mode
        toggleDisplayMode();
        return true;

      case K_ALT_R:  // repaint
        repaint();
        return true;

      case K_F6: {  // save all selected's
        FileSel fs;
        if (fs.Type("save selected fractals", "*").ExecuteSelectDir()) {
          setGeo();
          for (int r = 0; r < visGrid.GetCount(); r++) {
            auto _fg =
                ValueTo<FractalGeo>(visGrid.Get(r, 1));  // get center,range

            PNGEncoder().SaveFile(
                AppendFileName(fs.Get(), Format("frac-%d", r)),  // save
                genMandel(w, h, fg.iters, _fg.center, _fg.range));

            lsb =
                Format("saved file %d/%d", r, visGrid.GetCount());  // progress
            ProcessEvents();
          }
          lsb = Format("done, written %d files in %s ", visGrid.GetCount(),
                       fs.Get());
        }
      }
        return true;
      case K_F5: {  // save current fractal
        setGeo();

        FileSel fs;
        if (fs.Type("save PNG", "*.png").ExecuteSaveAs())
          PNGEncoder().SaveFile(fs.Get(),
                                genMandel(w, h, fg.iters, fg.center, fg.range));
      }
        return true;

      case K_ALT_N:  // new
        visGrid.Clear();

        fg.init();
        changed = true;
        repaint();
        return true;

      case K_ALT_LEFT:  // move
        fg.center -= Cmplx(deltaMove, 0.0);
        repaint();
        return true;
      case K_ALT_RIGHT:
        fg.center += Cmplx(deltaMove, 0.0);
        repaint();
        return true;
      case K_ALT_UP:
        fg.center -= Cmplx(0.0, deltaMove);
        repaint();
        return true;
      case K_ALT_DOWN:
        fg.center += Cmplx(0.0, deltaMove);
        repaint();
        return true;

      case K_ALT_ADD:  // + zoom in
        fg.range *= 0.8;
        repaint();
        return true;
      case K_ALT_SUBTRACT:  // zoom out
        fg.range /= 0.8;
        repaint();
        return true;

      case K_ENTER:  // redraw
        fg.iters = ~edIters;
        repaint();
        return true;
      case K_ALT_ENTER:  // save to selected
        changed = true;
        fg.iters = ~edIters;
        repaint();
        return true;
    }

    return TopWindow::Key(key, count);
  }

  void repaint() {
    int w = fracDisp.GetSize().cx, h = fracDisp.GetSize().cy;

    switch (dmode) {
      case Single: {
        imgFrac = genMandel(w, h, fg.iters, fg.center, fg.range);

        fracDisp.SetImage(imgFrac);

        if (changed) addImgGrid();

        lsb = Format(" lap:%0.f ms, size(%d,%d), |range|=%.1e", lap, w, h,
                     (double)abs(fg.range));
      } break;
      case Mosaic: {
      } break;
    }
    resizing = false;
  }

  void addImgGrid() {
    // add item to visGrid
    ImageCtrl *ic = new ImageCtrl;
    Image tn = genMandel(grSize, grSize, fg.iters, fg.center,
                         fg.range);  //  Rescale(imgFrac, Size(grSize, grSize));
    ic->SetImage(tn);

    int r = visGrid.GetCount();  // add fractal TN, center, range
    visGrid.Add("", RawToValue<FractalGeo>(fg), tn);

    visGrid.SetCtrl(r, 0, ic);
    visGrid.SetRowHeight(r + 1, grSize / 2);

    changed = false;
  }

  void reCalculateCenterRange(Point p) {
    Real w = fracDisp.GetSize().cx, h = fracDisp.GetSize().cy, dist = w / 2,
         rx = dist / w, ry = dist / h, ratio = abs(fg.range);

    fg.center += Cmplx(ratio * (w / 2 - p.x) / w, ratio * (h / 2 - p.y) / h);
    fg.range = Cmplx(fg.range.real() * rx, fg.range.imag() * ry);

    changed = true;
    repaint();
  }
};

GUI_APP_MAIN { FracWin().Run(); }
