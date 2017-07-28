#if defined(_MSC_VER)
#pragma once
#endif

#ifndef FS_TRIMESH_h
#define FS_TRIMESH_h

#include <vector>
#include "vec.h"
#include "types.h"

/************************************************************************/
/*				   Simple Triangle Mesh representation                  */
/************************************************************************/

struct TriMesh2d {
	std::vector<Vec3ui> tri;	// mesh connectivity
	std::vector<Vec2r> x;		// mesh vertex locations
	TriMesh2d() {
		clear();
	}

	~TriMesh2d() {
		tri.swap(vector<Vec3ui>(0));
		x.swap(vector<Vec2r>(0));
	}

	void clear() {
		tri.clear();
		x.clear();
	}

	void draw_mesh();

};

struct TriMesh {
	// list of triangles: the fundamental data
	std::vector<Vec3ui> tri;	// mesh connectivity
	std::vector<Vec3r> x;		// mesh vertex locations
	std::vector<Vec3r> normal;	// mesh normals

	std::vector<Real> masses;	// use for rigid body simulation,default to 0

	TriMesh() {
		tri.clear();
		x.clear();
		normal.clear();
	}

	~TriMesh() {
		clear();
	}

	void CopyTo(TriMesh &other) {
		other.tri.clear();
		other.x.clear();
		other.normal.clear();
		other.tri.resize(tri.size());
		other.x.resize(x.size());
		other.normal.resize(normal.size());
		other.tri = tri;
		other.x = x;
		other.normal = normal;
		//BLAS::copy(x, other.x);
		//BLAS::copy(tri, other.tri);
		//BLAS::copy(normal, other.normal);
	}
	
	void compute_normal(bool recompute = false);

	void draw_mesh(void) const;

	void scale(const Vec3r &center, Real length);

	void centerize();
	void scale(Real s);

	void flipXYZ();	// x->y y->z z->x
	void flipYZ();
	void flipXZ();
	void flipXY();

	// IO routines
	void dumpObj(const char *filename_format, ...);
	void dumpPbrt(const char *filename_format, ...);
	void dumpPovray(const char *filename_format, ...);

	void clear(void);
};


#endif // FS_TRIMESH_h
