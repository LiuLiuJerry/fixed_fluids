#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "util.h"
#include "trimesh.h"
#include "openglutils.h"

void TriMesh2d::draw_mesh() {
	draw_trimesh2d(x,tri);
}

void TriMesh::draw_mesh() const {
	if(normal.size()!=x.size())
		draw_trimesh3d(x,tri);
	else 
		draw_trimesh3d(x,tri,normal);
}

void TriMesh::clear(void) {
   tri.swap(vector<Vec3ui>(0));
   x.swap(vector<Vec3r>(0));
   normal.swap(vector<Vec3r>(0));
}

void TriMesh::compute_normal(bool recompute ) {
	if(normal.size() != 0 && !recompute)
		return;

	normal.resize(x.size(), Vec3r(0,0,0));
	
	for(unsigned int i = 0; i < tri.size(); i++) {
		int t0 = tri[i][0], t1 = tri[i][1], t2 = tri[i][2];
		Vec3r v0 = x[t0];
		Vec3r v1 = x[t1];
		Vec3r v2 = x[t2];
		Vec3r nml = cross(v1-v0,v2-v0);
		normal[t0] += nml;//*triangle_angle(v2,v0,v1);
		normal[t1] += nml;//*triangle_angle(v0,v1,v2);
		normal[t2] += nml;//*triangle_angle(v1,v2,v0);
	}
	// normalized them
	for(unsigned int i = 0; i < x.size(); i++) {
		if(mag(normal[i])>1e-12)
			normalize(normal[i]);
	}
}

void TriMesh::dumpObj( const char *filename_format, ... ) {
	char *filename=new char[256];
	va_list ap;
	va_start(ap, filename_format);

	vsprintf(filename, filename_format, ap);
	FILE *fp=fopen(filename, "wt");
	free(filename);
	va_end(ap);

	fprintf(fp, "#obj created by FluidSimulator written by Li Xiaosheng\n");
	
	// first vertex
	for(int i = 0; i < x.size(); i++) {
		fprintf(fp, "v %f %f %f\n", x[i][0], x[i][1], x[i][2]);
	}
	// face
	for(int i = 0; i < tri.size(); i++) {
		fprintf(fp, "f %d %d %d\n", tri[i][0]+1, tri[i][1]+1,tri[i][2]+1);
	}

	fclose(fp);
}

void TriMesh::dumpPbrt(const char *filename_format, ...) {
	char *filename=new char[256];
	va_list ap;
	va_start(ap, filename_format);

	vsprintf(filename, filename_format, ap);
	FILE *fp=fopen(filename, "wt");
	free(filename);
	va_end(ap);

	fprintf(fp, "# pbrt scene file created by FluidSimulator\n");
	fprintf(fp, "# @author : Li Xiaosheng\n");

	fprintf(fp,"AttributeBegin\n");
	fprintf(fp,"Shape \"trianglemesh\"\n");
	fprintf(fp,"\"integer indices\"\n");
	fprintf(fp,"[\n");
	for(int i = 0; i < tri.size(); i++)
		fprintf(fp, "	%d %d %d\n",tri[i][0],tri[i][1],tri[i][2]);
	fprintf(fp,"]\n");

	fprintf(fp,"\"point P\"\n");
	fprintf(fp,"[\n");
	for(int i = 0; i < x.size(); i++)
		fprintf(fp, "	%f %f %f\n",x[i][0],x[i][1],x[i][2]);
	fprintf(fp,"]\n");

	fprintf(fp,"\"normal N\"\n");
	fprintf(fp,"[\n");
	for(int i = 0; i < normal.size(); i++)
		fprintf(fp, "	%f %f %f\n",normal[i][0],normal[i][1],normal[i][2]);
	fprintf(fp,"]\n");

	fprintf(fp,"AttributeEnd\n");

	fclose(fp);
}

void TriMesh::dumpPovray(const char *filename_format, ...) {
	char *filename=new char[256];
	va_list ap;
	va_start(ap, filename_format);

	vsprintf(filename, filename_format, ap);
	FILE *fp=fopen(filename, "wt");
	free(filename);
	va_end(ap);

	fprintf(fp, "// pov-ray scene file created by FluidSimulator\n");
	fprintf(fp, "// @author : Liu Jiarui\n");
	fprintf(fp,"#include \"template.pov\"\n");
	fprintf(fp,"#include \"branch.pov\"\n");
	fprintf(fp,"#declare fluid_surface = mesh2 {\n");
	fprintf(fp,"vertex_vectors {\n");
	fprintf(fp,"%d,\n",x.size());
	int xsize = x.size();
	for(int i = 0; i < xsize; i++){
		fprintf(fp, "<%f, %f, %f>,\n",x[i][1],x[i][2],x[i][0]);
	}
	//fprintf(fp, "<%f, %f, %f>\n",x[x.size()-1][0],x[x.size()-1][1],x[x.size()-1][2]);
	fprintf(fp,"}\n");
	fprintf(fp,"normal_vectors {\n");
	fprintf(fp,"%d,\n", normal.size());
	for(int i = 0; i < normal.size(); i++){
		fprintf(fp, "<%f, %f, %f>,\n",normal[i][1], normal[i][2],normal[i][0]);
	}
	//fprintf(fp, "<%f, %f, %f>\n",normal[normal.size()-1][0],normal[normal.size()-1][1],normal[normal.size()-1][2]);
	fprintf(fp,"}\n");

	fprintf(fp,"face_indices {\n");
	fprintf(fp,"%d,\n", tri.size());
	for(int i = 0; i < tri.size(); i++){
		fprintf(fp, "<%d, %d, %d>,\n",tri[i][0],tri[i][1],tri[i][2]);
	}
	//fprintf(fp, "<%d, %d, %d>\n",tri[tri.size()-1][0],tri[tri.size()-1][1],tri[tri.size()-1][2]);
	fprintf(fp,"}\n");
	fprintf(fp,"}\n");
	fprintf(fp," #include \"footer.pov\" \n");

	fclose(fp);
}

void TriMesh::scale( const Vec3r &center, Real length ) {
	// find center (rough contain box)
	Vec3r min_corner(FLT_MAX,FLT_MAX,FLT_MAX), max_corner(-FLT_MAX,-FLT_MAX,-FLT_MAX);
	for(int i = 0; i < x.size(); i++) {
		min_corner[0] = min<Real>(x[i][0],min_corner[0]);
		min_corner[1] = min<Real>(x[i][1],min_corner[1]);
		min_corner[2] = min<Real>(x[i][2],min_corner[2]);

		max_corner[0] = max<Real>(x[i][0],max_corner[0]);
		max_corner[1] = max<Real>(x[i][1],max_corner[1]);
		max_corner[2] = max<Real>(x[i][2],max_corner[2]);
	}

	Vec3r old_center = (Real)0.5*(min_corner + max_corner);

	//Debug("old centre [%f,%f,%f]\n",old_center[0],old_center[1],old_center[2]);
	//Debug("new cenre [%f,%f,%f]\n",center[0],center[1],center[2]);

	Real maxLength = max<Real>(max_corner[0] - min_corner[0], max_corner[1] - min_corner[1], max_corner[2] - min_corner[2]);
	
	for(int i = 0; i < x.size(); i++) {
		x[i] = center + length*(x[i]-old_center)/maxLength;
	}
}

void TriMesh::scale(Real s) {
	for(int i = 0; i < x.size(); i++)
		x[i] *= s;
}

// force mass center to be origin
void TriMesh::centerize() {
	Vec3r center(0,0,0);
	for(int i = 0; i < x.size(); i++)
		center += x[i];
	center /= (Real)x.size();
	for(int i = 0; i < x.size(); i++)
		x[i] -= center;
}

void TriMesh::flipXYZ() {
	for(int i = 0; i < x.size(); i++) {
		Real temp = x[i][0];
		x[i][0] = x[i][1];
		x[i][1] = x[i][2];
		x[i][2] = temp;
		//swap(x[i][1],x[i][2]);
	}
}

void TriMesh::flipYZ() {
	for(int i = 0; i < x.size(); i++) {
		swap(x[i][1],x[i][2]);
	}
}

void TriMesh::flipXY() {
	for(int i = 0; i < x.size(); i++) {
		swap(x[i][1],x[i][0]);
	}
}

void TriMesh::flipXZ() {
	for(int i = 0; i < x.size(); i++) {
		swap(x[i][0],x[i][2]);
	}
}