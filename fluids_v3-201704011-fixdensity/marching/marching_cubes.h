#if defined(_MSC_VER)
#pragma once
#endif

#ifndef FS_MARCHING_CUBES_h
#define FS_MARCHING_CUBES_h

/**
 * @file    MarchingCubes.h
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.2
 * @date    12/08/2002
 *
 * @brief   MarchingCubes Algorithm, modified by Li Xiaosheng - 2011-9-18
 */
//________________________________________________

#include "types.h"
#include "trimesh.h"
#include "array3.h"


//_____________________________________________________________________________
/** Marching Cubes algorithm wrapper */
/** \class MarchingCubes
  * \brief Marching Cubes algorithm.
  */
class MarchingCubes
//-----------------------------------------------------------------------------
{
// Constructors
public :
  /**
   * Main and default constructor
   * \brief constructor
   * \param size_x width  of the grid
   * \param size_y depth  of the grid
   * \param size_z height of the grid
   */
  MarchingCubes ( const Vec3r &origin, real dx, const Array3r &phi);
  /** Destructor */
  ~MarchingCubes() ;

//-----------------------------------------------------------------------------
// Accessors
public :
  /** accesses the number of vertices of the generated mesh */
  inline const int nverts() const { return trimesh.x.size() ; }
  /** accesses the number of triangles of the generated mesh */
  inline const int ntrigs() const { return trimesh.tri.size() ; }

  /**  accesses the width  of the grid */
  inline const int size_x() const { return _size_x ; }
  /**  accesses the depth  of the grid */
  inline const int size_y() const { return _size_y ; }
  /**  accesses the height of the grid */
  inline const int size_z() const { return _size_z ; }
  /**
   * selects wether the algorithm will use the enhanced topologically controlled lookup table or the original MarchingCubes
   * \param originalMC true for the original Marching Cubes
   */
  inline void set_method    ( const bool originalMC = false ) { _originalMC = originalMC ; }

  // Data access
  /**
   * accesses a specific cube of the grid
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline const real get_data  ( const int i, const int j, const int k ) const { return phi(i,j,k); }
  /**
   * sets a specific cube of the grid
   * \param val new value for the cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
 
  // Data initialization
  /** inits all structures (must set sizes before call) : the temporary structures and the mesh buffers */
  void init_all   () ;

  /** clears all structures : the temporary structures and the mesh buffers */
  void clean_all  () ;

  void draw_mesh();

//-----------------------------------------------------------------------------
// Algorithm
public :
  /**
   * Main algorithm : must be called after init_all
   * \param iso isovalue
   */

  void run( real iso = (real)0.0 ) ;
  void contour_grid(real iso = (real)0.0) { 
	  init_all();
	  run(iso);
  }
  void swap_trimesh(TriMesh &_trimesh);
  inline TriMesh* get_mesh() {
	  return &trimesh;
  }

protected :
  /** tesselates one cube */
  void process_cube ()             ;
  /** tests if the components of the tesselation of the cube should be connected by the interior of an ambiguous face */
  bool test_face    ( schar face ) ;
  /** tests if the components of the tesselation of the cube should be connected through the interior of the cube */
  bool test_interior( schar s )    ;


//-----------------------------------------------------------------------------
// Operations
protected :
  /**
   * computes almost all the vertices of the mesh by interpolation along the cubes edges
   * \param iso isovalue
   */
  void compute_intersection_points( real iso ) ;

  /**
   * routine to add a triangle to the mesh
   * \param trig the code for the triangle as a sequence of edges index
   * \param n    the number of triangles to produce
   * \param v12  the index of the interior vertex to use, if necessary
   */
  void add_triangle ( const char* trig, char n, int v12 = -1 ) ;

  /** tests and eventually doubles the vertex buffer capacity for a new vertex insertion */
  void test_vertex_addition() ;
  /** adds a vertex on the current horizontal edge */
  int add_x_vertex() ;
  /** adds a vertex on the current longitudinal edge */
  int add_y_vertex() ;
  /** adds a vertex on the current vertical edge */
  int add_z_vertex() ;
  /** adds a vertex inside the current cube */
  int add_c_vertex() ;

  /**
   * interpolates the horizontal gradient of the implicit function at the lower vertex of the specified cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  real get_x_grad( const int i, const int j, const int k ) const ;
  /**
   * interpolates the longitudinal gradient of the implicit function at the lower vertex of the specified cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  real get_y_grad( const int i, const int j, const int k ) const ;
  /**
   * interpolates the vertical gradient of the implicit function at the lower vertex of the specified cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  real get_z_grad( const int i, const int j, const int k ) const ;

  /**
   * accesses the pre-computed vertex index on the lower horizontal edge of a specific cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline int   get_x_vert( const int i, const int j, const int k ) const { return _x_verts[ i + j*_size_x + k*_size_x*_size_y] ; }
  /**
   * accesses the pre-computed vertex index on the lower longitudinal edge of a specific cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline int   get_y_vert( const int i, const int j, const int k ) const { return _y_verts[ i + j*_size_x + k*_size_x*_size_y] ; }
  /**
   * accesses the pre-computed vertex index on the lower vertical edge of a specific cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline int   get_z_vert( const int i, const int j, const int k ) const { return _z_verts[ i + j*_size_x + k*_size_x*_size_y] ; }

  /**
   * sets the pre-computed vertex index on the lower horizontal edge of a specific cube
   * \param val the index of the new vertex
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_x_vert( const int val, const int i, const int j, const int k ) { _x_verts[ i + j*_size_x + k*_size_x*_size_y] = val ; }
  /**
   * sets the pre-computed vertex index on the lower longitudinal edge of a specific cube
   * \param val the index of the new vertex
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_y_vert( const int val, const int i, const int j, const int k ) { _y_verts[ i + j*_size_x + k*_size_x*_size_y] = val ; }
  /**
   * sets the pre-computed vertex index on the lower vertical edge of a specific cube
   * \param val the index of the new vertex
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_z_vert( const int val, const int i, const int j, const int k ) { _z_verts[ i + j*_size_x + k*_size_x*_size_y] = val ; }

  /** prints cube for debug */
  void print_cube() ;

//-----------------------------------------------------------------------------
// Elements
protected :
  bool      _originalMC ;   /**< selects wether the algorithm will use the enhanced topologically controlled lookup table or the original MarchingCubes */
  bool      _ext_data   ;   /**< selects wether to allocate data or use data from another class */

  int       _size_x     ;  /**< width  of the grid */
  int       _size_y     ;  /**< depth  of the grid */
  int       _size_z     ;  /**< height of the grid */
  const Array3r	&phi	;  /**< implicit function values sampled on the grid */ 
  Vec3r		 origin		;
  real		 dx			;

  int      *_x_verts    ;  /**< pre-computed vertex indices on the lower horizontal   edge of each cube */
  int      *_y_verts    ;  /**< pre-computed vertex indices on the lower longitudinal edge of each cube */
  int      *_z_verts    ;  /**< pre-computed vertex indices on the lower vertical     edge of each cube */

  TriMesh	trimesh		;

  int       _i          ;  /**< abscisse of the active cube */
  int       _j          ;  /**< height of the active cube */
  int       _k          ;  /**< ordinate of the active cube */

  real      _cube[8]    ;  /**< values of the implicit function on the active cube */
  uchar     _lut_entry  ;  /**< cube sign representation in [0..255] */
  uchar     _case       ;  /**< case of the active cube in [0..15] */
  uchar     _config     ;  /**< configuration of the active cube */
  uchar     _subconfig  ;  /**< subconfiguration of the active cube */
};
//_____________________________________________________________________________



#endif // FS_MARCHING_CUBES_h



