#include <CtrlLib/CtrlLib.h>
#include <GLCtrl/GLCtrl.h>
#include "poly/parser.hpp"

using namespace Upp;

class PolyWin : public TopWindow {
  struct GLPoly : GLCtrl {
    Point point;
    float zoom = 2;
    Polyhedron poly;

    float anglex = 0, angley = 0;
    Point last;
    bool down = false;

    void setPoly(Polyhedron &poly) {
      this->poly = poly;
      Refresh();
    }

    void init() {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glLoadIdentity();
    }

    void setRotate() {
      glTranslatef(0, 0, -zoom);

      glRotatef(anglex, 1, 0, 0);
      glRotatef(-angley, 0, 1, 0);
    }

    void GLPaint() override {
      StdView();

      init();
      setRotate();

      int f = 0;  // current face index
      for (auto face : poly.faces) {
        glBegin(GL_POLYGON);
        glColor3fv((GLfloat *)&poly.colors[f]);

        for (auto ix : face) {
          glNormal3fv((GLfloat *)&poly.normals[ix]);
          glVertex3fv((GLfloat *)&poly.vertexes[ix]);
        }
        f++;
        glEnd();
      }

      if (poly.faces.size() < 500) {
        glColor3f(1, 0, 0);
        for (auto face : poly.faces) {  // draw line
          glBegin(GL_LINE_LOOP);
          for (auto ix : face) glVertex3fv((GLfloat *)&poly.vertexes[ix]);
          glEnd();
        }
      }
    }

    void LeftDown(Point p, dword) override { down = true; }
    void LeftUp(Point p, dword) override { down = false; }

    void RightDown(Point p, dword) override {
      poly.calc_colors();
      Refresh();
    }

    void MouseMove(Point p, dword flags) override {
      if (down) {
        float d = 2;

        angley += (p.x - last.x) / d;
        anglex += (p.y - last.y) / d;

        last = p;
        Refresh();
      }
    }
    void MouseWheel(Point p, int zdelta, dword) override {
      zoom += zdelta / 180.0;
      Refresh();
    }
  };

  GLPoly glPoly;
  TimeStop ts;

  StatusBar sb;
  EditString edPoly;
  Label lsb;

  Polyhedron poly;
  double lap;

 public:
  PolyWin() {
    Title("polyhedronisme").Sizeable().SetRect(Size(2000, 1500));

    Add(glPoly.SizePos());

    edPoly.SetText(Parser::randomTransform());
    edPoly.Tip("transformations: dagprPqkcwnxlH, poly: TCIODPAYUVJ#");
    edPoly.SetRect(0, 0, 210, 50);
    sb.Add(edPoly);

    lsb.SetRect(220, 0, 2000, 50);
    sb.Add(lsb);

    Key(K_ENTER, 0);

    AddFrame(sb);
  }

  void msg() {
    lsb = Format("%d faces, %d vertex, lap:%.0f ms", (int)poly.n_faces(),
                 (int)poly.n_vertexes(), lap / 1000);
  }

  bool Key(dword key, int flags) override {
    switch (key) {
      case K_ESCAPE:
        Close();

      case K_ENTER:
        ts.Reset();

        poly = Parser::parse((~edPoly).ToStd());
        lap = ts.Elapsed();

        glPoly.setPoly(poly);

        msg();

        return true;

      case K_F9:
      case K_ALT_R:
        ts.Reset();

        poly = Parser::randomPoly();
        glPoly.setPoly(poly);

        edPoly.SetText(poly.name);
        lap = ts.Elapsed();

        msg();

        return true;
    }
    return TopWindow::Key(key, flags);
  }
};

GUI_APP_MAIN {
  srand(time(nullptr));

  // Ctrl::GlobalBackPaint();

  PolyWin().Run();
}

