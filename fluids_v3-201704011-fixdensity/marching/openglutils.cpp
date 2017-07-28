#include "openglutils.h"

#ifdef __APPLE__
#include <GLUT/glut.h> // why does Apple have to put glut.h here...
#else
#include <GL/glut.h> // ...when everyone else puts it here?
#endif

#include "vec.h"
#include <cfloat>

void draw_circle2d(const Vec2r& centre, Real rad, int segs)
{
	glLineWidth(2);
   glBegin(GL_POLYGON);
   for(int i=0;i<segs;i++){
	  Real cosine=rad*cos(i*2*M_PI/(Real)(segs));
	  Real sine=rad* sin(i*2*M_PI/(Real)(segs));
	  glVertex2f((float)(cosine+centre[0]),(float)(sine+centre[1]));
	 // glVertex2fv((Vec2f(cosine,sine) + centre).v);
   }
   glEnd();
}
void draw_circle(const Vec2r& centre, Real rad, int segs)
{
   //glLineWidth(2);
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   glBegin(GL_POLYGON);
   //glColor3f(0,1,0);
   for(int i=0;i<segs;i++){
	  Real cosine=rad*cos(i*2*M_PI/(Real)(segs));
	  Real sine=rad* sin(i*2*M_PI/(Real)(segs));
	  glVertex2f((float)(cosine+centre[0]),(float)(sine+centre[1]));
	 // glVertex2fv((Vec2f(cosine,sine) + centre).v);
   }
   glEnd();
}

void draw_grid2d(const Vec2r& origin, Real dx, int nx, int ny) {
   Real width = nx*dx;
   Real height = ny*dx;
   glColor3f(0.753, 0, 0.753);
   glLineWidth(0.5);
   glBegin(GL_LINES);
   for(int i = 0; i <= nx; i++) {
	  Vec2r a(i*dx, 0);
	  Vec2r b(i*dx, height);
	  glVertex2f(origin[0]+a[0],origin[1]+a[1]);
	  glVertex2f(origin[0]+b[0],origin[1]+b[1]);
	 // glVertex2fv((origin+a).v); 
	 // glVertex2fv((origin+b).v);
   }
   for(int j = 0; j <= ny; ++j) {
	  Vec2r a(0,j*dx);
	  Vec2r b(width,j*dx);
	  glVertex2f(origin[0]+a[0],origin[1]+a[1]);
	  glVertex2f(origin[0]+b[0],origin[1]+b[1]);
	 // glVertex2fv((origin + a).v); 
	  //glVertex2fv((origin + b).v);
   }
   glEnd();
}

void draw_box2d(const Vec2r& origin, Real width, Real height) {
   glBegin(GL_POLYGON);
   /*glVertex2fv(origin.v);
   glVertex2fv((origin + Vec2r(0, height)).v);
   glVertex2fv((origin + Vec2r(width, height)).v);
   glVertex2fv((origin + Vec2r(width, 0)).v);*/
   glVertex2f(origin[0],origin[1]);
   glVertex2f(origin[0],origin[1]+height);
   glVertex2f(origin[0]+width,origin[1]+height);
   glVertex2f(origin[0]+width,origin[1]);
   glEnd();
}

void draw_segmentset2d(const std::vector<Vec2r>& vertices, const std::vector<Vec2ui>& edges) {
   glBegin(GL_LINES);
   for(unsigned int i = 0; i < edges.size(); ++i) {
	 /* glVertex2fv(vertices[edges[i][0]].v);      
	  glVertex2fv(vertices[edges[i][1]].v);*/
	   glVertex2f(vertices[edges[i][0]][0],vertices[edges[i][0]][1]);
	   glVertex2f(vertices[edges[i][1]][0],vertices[edges[i][1]][1]);
   }
   glEnd();
}

void draw_segmentset2d(const std::vector<Vec2r>& vertices, const std::vector<Vec2i>& edges) {
   glBegin(GL_LINES);
   for(unsigned int i = 0; i < edges.size(); ++i) {
	 /* glVertex2fv(vertices[edges[i][0]].v);      
	  glVertex2fv(vertices[edges[i][1]].v);*/
	   glVertex2f(vertices[edges[i][0]][0],vertices[edges[i][0]][1]);
	   glVertex2f(vertices[edges[i][1]][0],vertices[edges[i][1]][1]);
   }
   glEnd();
}

void draw_points2d(const std::vector<Vec2r>& points) {
   glBegin(GL_POINTS);
   for(unsigned int i = 0; i < points.size(); ++i) {
	  glVertex2f(points[i][0],points[i][1]);      
   }
   glEnd();
}
void draw_points2d_color(const std::vector<Vec2r>& points,const std::vector<Vec3r>& color,Real radius) {
  // glBegin(GL_POINTS);
   for(unsigned int i = 0; i < points.size(); ++i) {
	  glColor3f(color[i][0],color[i][1],color[i][2]);
	  //glVertex2fv(points[i].v);
	  draw_circle(points[i],radius,10);
   }
 //  glEnd();
}
void draw_points_cloud(const std::vector<Vec2r>& points,Real radius) {
  // glBegin(GL_POINTS);
   for(unsigned int i = 0; i < points.size(); ++i) {
	 // glColor3fv(color[i].v);
	 // glVertex2fv(points[i].v);
	   draw_circle(points[i],radius,10);	   
   }
 //  glEnd();
}
void draw_polygon2d(const std::vector<Vec2r>& vertices) {
   glBegin(GL_POLYGON);
   for(unsigned int i = 0; i < vertices.size(); ++i)
	  glVertex2f(vertices[i][0],vertices[i][1]);      
   glEnd();
}

void draw_polygon2d(const std::vector<Vec2r>& vertices, const std::vector<int>& order) {
   glBegin(GL_POLYGON);
   for(unsigned int i = 0; i < order.size(); ++i)
	  glVertex2f(vertices[order[i]][0],vertices[order[i]][1]);      
   glEnd();

}
void draw_segment2d(const Vec2r& start, const Vec2r& end) {
   glBegin(GL_LINES);
   glVertex2f(start[0],start[1]);      
   glVertex2f(end[0],end[1]);      
   glEnd();
}

void draw_arrow2d(const Vec2r& start, const Vec2r& end, Real arrow_head_len)
{
   Vec2r direction = end - start;

   Vec2r dir_norm = direction;
   
   //TODO Possibly automatically scale arrowhead length based on vector magnitude
   if(mag(dir_norm) < 1e-14)
	  return;

   normalize(dir_norm);
   Vec2r perp(dir_norm[1],-dir_norm[0]);

   Vec2r tip_left = end + arrow_head_len/(Real)sqrt(2.0)*(-dir_norm + perp);
   Vec2r tip_right = end + arrow_head_len/(Real)sqrt(2.0)*(-dir_norm - perp);

   glColor3f(1,1,0);
   glLineWidth(1);
   
   glBegin(GL_LINES);
   glVertex2f(start[0],start[1]);
   glVertex2f(end[0],end[1]);
   glVertex2f(end[0],end[1]);
   glVertex2f(tip_left[0],tip_left[1]);
   glVertex2f(end[0],end[1]);
   glVertex2f(tip_right[0],tip_right[1]);
   glEnd();
	
}

void draw_arrow3d(const Vec3r& start, const Vec3r& end, Real arrow_head_len) {
	Vec3r direction = end - start;

	Vec3r dir_norm = direction;

	//TODO Possibly automatically scale arrowhead length based on vector magnitude
	if(mag(dir_norm) < 1e-14)
		return;

	normalize(dir_norm);
	Vec3r perp(1,-dir_norm[0]/dir_norm[1],0);

	Vec3r tip_left = end + arrow_head_len/(Real)sqrt(2.0)*(-dir_norm + perp);
	Vec3r tip_right = end + arrow_head_len/(Real)sqrt(2.0)*(-dir_norm - perp);

	glColor3f(1,1,0);
	glLineWidth(1);

	glBegin(GL_LINES);
	glVertex3f(start[0],start[1],start[2]);
	glVertex3f(end[0],end[1],end[2]);
	glVertex3f(end[0],end[1],end[2]);
	glVertex3f(tip_left[0],tip_left[1],tip_left[2]);
	glVertex3f(end[0],end[1],end[2]);
	glVertex3f(tip_right[0],tip_right[1],tip_right[2]);
	glEnd();
}

void draw_trimesh2d(const std::vector<Vec2r>& vertices, const std::vector<Vec3ui>& tris) {
   glBegin(GL_TRIANGLES);
   for(unsigned int i = 0; i < tris.size(); ++i) {
	  glVertex2f(vertices[tris[i][0]][0],vertices[tris[i][0]][1]);
	  glVertex2f(vertices[tris[i][1]][0],vertices[tris[i][1]][1]);
	  glVertex2f(vertices[tris[i][2]][0],vertices[tris[i][2]][1]);
	}
   glEnd();
}
	   
	   
void hueToRGB(Real hue, Real sat, Real val, Real &r, Real &g, Real &b) {   
   //compute hue (adapted from an older Wikipedia article)
   int Hi = (int)(floor(hue / 60.0f)) % 6;
   Real f = hue / 60 - Hi;
   Real p = val * (1 - sat);
   Real q = val * (1- f * sat);
   Real t = val * (1 - (1 - f) * sat);
   
   switch(Hi) {
	  case 0:
		 r=val;
		 g=t;
		 b=p;
		 break;
	  case 1:
		 r=q;
		 g=val;
		 b=p;
		 break;
	  case 2:
		 r=p;
		 g=val;
		 b=t;
		 break;
	  case 3:
		 r=p;
		 g=q;
		 b=val;
		 break;
	  case 4:
		 r=t;
		 g=p;
		 b=val;
		break;
	  case 5:
		 r=val;
		 g=p;
		 b=q;
		 break;
   }
}

void draw_grid_data2d(Array2r& data, Vec2r origin, Real dx, bool color) {
   Real max_val = FLT_MIN;
   Real min_val = FLT_MAX;
   for(int j = 0; j < data.nj; ++j) for(int i = 0; i < data.ni; ++i) {
	  max_val = max(data(i,j),max_val);
	  min_val = min(data(i,j),min_val);
   }
   
   for(int j = 0; j < data.nj; ++j) {
	  for(int i = 0; i < data.ni; ++i) {
		 glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		 Vec2r bl = origin + Vec2r(i*dx,j*dx);
		 Real r,g,b;
		 if(color) {
			hueToRGB(240*(data(i,j) - min_val)/(max_val-min_val), 1, 1, r,g,b);
		 }
		 else {
			Real gray = (data(i,j) - min_val)/(max_val-min_val);
			r = g = b = gray;
		 }
		 //TODO Black body colormap, if I can find it.
		 glColor3f(r,g,b);
		 draw_box2d(bl, dx, dx);
	  }
   }

}

void draw_trimesh3d(const std::vector<Vec3r>& vertices, const std::vector<Vec3ui>& tris) {
   
	std::vector<Vec3r> normals;
	normals.resize(vertices.size());
	for(unsigned int i = 0; i < tris.size(); i++) {
		Vec3r v1 = vertices[tris[i][1]] - vertices[tris[i][0]]; 
		Vec3r v2 = vertices[tris[i][2]] - vertices[tris[i][0]];
		Vec3r normal = Vec3r(v1[1]*v2[2]-v1[2]*v2[1],v2[0]*v1[2]-v1[0]*v2[2],v1[0]*v2[1]-v2[0]*v1[1]);
		normalize(normal);
		normals[tris[i][0]] += normal;
		normals[tris[i][1]] += normal;
		normals[tris[i][2]] += normal;
	}

	for(unsigned int i = 0; i < normals.size(); i++)
		normalize(normals[i]);



	draw_trimesh3d(vertices, tris, normals);

#if 0
   glBegin(GL_TRIANGLES);
   for(unsigned int i = 0; i < tris.size(); ++i) {

	 /*  Vec3r v1 = vertices[tris[i][1]] - vertices[tris[i][0]]; 
	   Vec3r v2 = vertices[tris[i][2]] - vertices[tris[i][0]];
	   Vec3r normal = Vec3r(v1[1]*v2[2]-v1[2]*v2[1],v2[0]*v1[2]-v1[0]*v2[2],v1[0]*v2[1]-v2[0]*v1[1]);
	   normalize(normal);
	   glNormal3f(normal[0],normal[1],normal[2]);*/
	  glVertex3f(vertices[tris[i][0]][0],vertices[tris[i][0]][1],vertices[tris[i][0]][2]);
	  glVertex3f(vertices[tris[i][1]][0],vertices[tris[i][1]][1],vertices[tris[i][1]][2]);
	  glVertex3f(vertices[tris[i][2]][0],vertices[tris[i][2]][1],vertices[tris[i][2]][2]);
   }
   glEnd();
#endif
}

void draw_trimesh3d(const std::vector<Vec3r>& vertices, const std::vector<Vec3ui>& tris, const std::vector<Vec3r> & normals) {
   glBegin(GL_TRIANGLES);
  // glNormal3f(0,0,1);
   for(unsigned int i = 0; i < tris.size(); ++i) {
	  glNormal3f(normals[tris[i][0]][0],normals[tris[i][0]][1],normals[tris[i][0]][2]);
	  glVertex3f(vertices[tris[i][0]][0],vertices[tris[i][0]][1],vertices[tris[i][0]][2]);
	  glNormal3f(normals[tris[i][1]][0],normals[tris[i][1]][1],normals[tris[i][1]][2]);
	  glVertex3f(vertices[tris[i][1]][0],vertices[tris[i][1]][1],vertices[tris[i][1]][2]);
	  glNormal3f(normals[tris[i][2]][0],normals[tris[i][2]][1],normals[tris[i][2]][2]);
	  glVertex3f(vertices[tris[i][2]][0],vertices[tris[i][2]][1],vertices[tris[i][2]][2]);
   }
   glEnd();
}

void draw_box3d(const Vec3r& dimensions) {
   
   //Draw an axis-aligned box with specified dimensions, 
   //where the midpoint of the box is at the origin

   Real width = dimensions[0];
   Real height = dimensions[1];
   Real depth = dimensions[2];

   glPushMatrix();
   glTranslatef(width, height, depth);

   glBegin(GL_POLYGON);
   glNormal3f(-1,0,0);
   glVertex3f(-0.5*width, -0.5*height, 0.5*depth);
   glVertex3f(-0.5*width, 0.5*height, 0.5*depth);
   glVertex3f(-0.5*width, 0.5*height, -0.5*depth);
   glVertex3f(-0.5*width, -0.5*height, -0.5*depth);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(1,0,0);
   glVertex3f(0.5*width, -0.5*height, 0.5*depth);
   glVertex3f(0.5*width, 0.5*height, 0.5*depth);
   glVertex3f(0.5*width, 0.5*height, -0.5*depth);
   glVertex3f(0.5*width, -0.5*height, -0.5*depth);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(0,0,-1);
   glVertex3f(-0.5*width, -0.5*height, -0.5*depth);
   glVertex3f(0.5*width, -0.5*height, -0.5*depth);
   glVertex3f(0.5*width, 0.5*height, -0.5*depth);
   glVertex3f(-0.5*width, 0.5*height, -0.5*depth);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(0,0,1);
   glVertex3f(-0.5*width, -0.5*height, 0.5*depth);
   glVertex3f(0.5*width, -0.5*height, 0.5*depth);
   glVertex3f(0.5*width, 0.5*height, 0.5*depth);
   glVertex3f(-0.5*width, 0.5*height, 0.5*depth);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(0,-1,0);
   glVertex3f(-0.5*width, -0.5*height, 0.5*depth);
   glVertex3f(0.5*width, -0.5*height, 0.5*depth);
   glVertex3f(0.5*width, -0.5*height, -0.5*depth);
   glVertex3f(-0.5*width, -0.5*height, -0.5*depth);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(0,1,0);
   glVertex3f(-0.5*width, 0.5*height, 0.5*depth);
   glVertex3f(0.5*width, 0.5*height, 0.5*depth);
   glVertex3f(0.5*width, 0.5*height, -0.5*depth);
   glVertex3f(-0.5*width, 0.5*height, -0.5*depth);
   glEnd();
   glPopMatrix();
}

void draw_line_box3d(const Vec3r& dimensions)
{
	Real width = dimensions[0];
	Real height = dimensions[1];
	Real depth = dimensions[2];

	glDisable(GL_LIGHTING);
	glColor3f(0,1,1);
	glBegin(GL_LINES);

	// front
	glVertex3f(0,0,0); glVertex3f(width,0,0);
	glVertex3f(0,0,0); glVertex3f(0,height,0);
	glVertex3f(width,height,0); glVertex3f(width,0,0);
	glVertex3f(width,height,0); glVertex3f(0,height,0);
	// back
	glVertex3f(0,0,depth); glVertex3f(width,0,depth);
	glVertex3f(0,0,depth); glVertex3f(0,height,depth);
	glVertex3f(width,height,depth); glVertex3f(width,0,depth);
	glVertex3f(width,height,depth); glVertex3f(0,height,depth);
	// left
	glVertex3f(0,0,0); glVertex3f(0,0,depth);
	glVertex3f(0,height,0); glVertex3f(0,height,depth);
	// right
	glVertex3f(width,0,0); glVertex3f(width,0,depth);
	glVertex3f(width,height,0); glVertex3f(width,height,depth);

	glEnd();

	glEnable(GL_LIGHTING);
}

void draw_line_box3d(const Vec3r& left, const Vec3r& right) {

}

void draw_densities(const Vec2r& position, Real dx,Real p1,Real p2,Real p3,Real p4)
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_POLYGON);
	glColor3f(p1,p1,p1);glVertex2f(position[0],position[1]);
	glColor3f(p2,p2,p2);glVertex2f(position[0]+dx,position[1]);
	glColor3f(p3,p3,p3);glVertex2f(position[0]+dx,position[1]+dx);
	glColor3f(p4,p4,p4);glVertex2f(position[0],position[1]+dx);
	glEnd();
}

void drawBitmapString( const char *string) {
	while (*string) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *string++);
}
