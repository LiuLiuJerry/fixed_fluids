//---------------------------------------------------------------------------//
// AccessObj.h: interface for the CAccessObj class.                          //
//---------------------------------------------------------------------------//
#ifndef __wtypes_h__
#include <wtypes.h>
#endif

#ifndef __WINDEF_
#include <windef.h>
#endif
//#include "assert.h"
#include "../common/app_util.h"

#ifndef _MY_ACCESS_OBJ_
#define _MY_ACCESS_OBJ_
//#define CellSize 10

#define OBJ_NONE (0)		                       // render with only vertices
#define OBJ_FLAT (1 << 0)	                       // render with facet normals
#define OBJ_SMOOTH (1 << 1)	                      // render with vertex normals
#define OBJ_TEXTURE (1 << 2)	                  // render with texture coords
#define OBJ_COLOR (1 << 3)	                              // render with colors
#define OBJ_MATERIAL (1 << 4)                          // render with materials

#define objMax(a,b)	(((a)>(b))?(a):(b))
#define objMin(a,b)	(((a)<(b))?(a):(b))
#define objAbs(x)	(((x)>0.f)?(x):(-x))

#define T(x) (m_pModel->pTriangles[(x)])

struct node{
	int TriNum;//一个node代表一个三角形的记录，TriNum是这个三角形的编号
	struct node * next;
};

//---------------------------------------------------------------------------//
// Definition of the OBJ R/W class                                           //
//---------------------------------------------------------------------------//
class CAccessObj
{
public:
	// --------------------------------------------------------------------
	// COBJmaterial: defines a material in a model. 
	class COBJmaterial
	{
	public:
		char name[256];			// name of material 
		float diffuse[4];		// diffuse component 
		float ambient[4];		// ambient component 
		float specular[4];		// specular component 
		float emissive[4];		// emissive component 
		float shininess[1];		// specular exponent 
		
		char * MaterialName() {return (name);}
		float * DiffuseComponent() {return (diffuse);}
		float * AmbientComponent() {return (ambient);}
		float * SpecularComponent() {return (specular);}
		float * EmissiveComponent() {return (emissive);}
		float * ShininessComponent() {return (shininess);}
		
		COBJmaterial()
		{
			sprintf (name, "default");;
			diffuse[0] = diffuse[1] = diffuse[2] = diffuse[3] = 1.0f;
			ambient[0] = ambient[1] = ambient[2] = ambient[3] = 0.1f;
			specular[0] = specular[1] = specular[2] = specular[3] = 0.0f;
			emissive[0] = emissive[1] = emissive[2] = emissive[3] = 0.0f;
			shininess[0] = 10;
		}
		
		virtual ~COBJmaterial() {}
	};

	// ----------------------------------------------------------------------//
	// COBJtriangle: defines a triangle in a model.                          //
	// ----------------------------------------------------------------------//

	class COBJtriangle
	{
	public:
		unsigned int VertIdx[3];    // array of triangle vertex indices 
		unsigned int NormIdx[3];    // array of triangle normal indices 
		unsigned int TexIdx[3];     // array of triangle texcoord indices
		unsigned int findex;		// index of triangle facet normal 
		COBJtriangle()	{}
		virtual ~COBJtriangle()	{}
	};

	// ----------------------------------------------------------------------//
	// COBJgroup: defines a group in a model.                                //
	// ----------------------------------------------------------------------//

	class COBJgroup
	{
	public:
		char name[256];			                          // name of this group
		unsigned int nTriangles;	                      // number of triangles in this group
		unsigned int* pTriangles;	                      // array of triangle indices
		unsigned int material;		                      // index to material for group
		int material_id;
		char texturename[256];
		bool m_bTexture;
		int m_iTextureType;
		float m_fTran_X;
		float m_fTran_Y;
		float m_fTran_Z;
		float m_fScale_X;
		float m_fScale_Y;
		float m_fScale_Z;
		float m_fRotate_X;
		float m_fRotate_Y;
		float m_fRotate_Z;
		class COBJgroup* next;		                       // pointer to next group in model
		
		COBJgroup()
		{
			nTriangles = 0;
			pTriangles = NULL;
			next = NULL;
		}
		
		virtual ~COBJgroup()
		{
			if (nTriangles != 0)
			{
				delete [] pTriangles;
				pTriangles = NULL;
			}
		}
	};

	// ----------------------------------------------------------------------//
	// COBJmodel: defines a model.                                           //
	// ----------------------------------------------------------------------//
	class COBJmodel
	{
	public:
		char pathname[256];	        // path to this model 
		char mtllibname[256];		// name of the material library 
		unsigned int nVertices;		// number of vertices in model 
		Vector3DF* vpVertices;		// array of vertices 
		unsigned int nNormals;		// number of normals in model 
		Vector3DF* vpNormals;		    // array of normals 
		unsigned int nTexCoords;	// number of texcoords in model 
		Vector3DF* vpTexCoords;		// array of texture coordinates 
		unsigned int nFacetnorms;	// number of facetnorms in model 
		Vector3DF* vpFacetNorms;		// array of facetnorms 
		unsigned int nTriangles;	// number of triangles in model 
		COBJtriangle* pTriangles;	// array of triangles 
		unsigned int nMaterials;	// number of materials in model 
		COBJmaterial* pMaterials;	// array of materials 
		unsigned int nGroups;		// number of groups in model 
		COBJgroup* pGroups;		    // linked list of groups 
		Vector3DF position;		    // position of the model 
		int  m_nMaterial;


		// construction
		COBJmodel()
		{
			nVertices   = 0;
			vpVertices  = NULL;
			nNormals    = 0;
			vpNormals   = NULL;
			nTexCoords  = 0;
			vpTexCoords = NULL;
			nFacetnorms = 0;
			vpFacetNorms= NULL;
			nTriangles  = 0;
			pTriangles  = NULL;
			nMaterials  = 0;
			pMaterials  = NULL;
			nGroups     = 0;
			pGroups     = NULL;
			position    = Vector3DF(0, 0, 0);
			m_nMaterial = 0;

		}

		// free all memory
		void	Destroy()
		{
			COBJgroup *group;
			
			if (vpVertices)		delete [] vpVertices;
			if (vpNormals)		delete [] vpNormals;
			if (vpTexCoords)	delete [] vpTexCoords;
			if (vpFacetNorms)	delete [] vpFacetNorms;
			if (pTriangles)
			{
				delete [] pTriangles;
			}
			if (pMaterials)		delete [] pMaterials;
			
			while(pGroups)
			{
				group = pGroups;
				pGroups = pGroups->next;
				delete group;
			}
			
			nVertices    = 0;
			vpVertices   = NULL;
			nNormals     = 0;
			vpNormals    = NULL;
			nTexCoords   = 0;
			vpTexCoords  = NULL;
			nFacetnorms  = 0;
			vpFacetNorms = NULL;
			nTriangles   = 0;
			pTriangles   = NULL;
			nMaterials   = 0;
			pMaterials   = NULL;
			nGroups      = 0;
			pGroups      = NULL;
			position     = Vector3DF(0, 0, 0);
		}
		
		// destruction
		virtual ~COBJmodel()
		{
			Destroy();
		}
	};

	//-----------------------------------------------------------------------//
	// A temporal calss                                                      //
	//-----------------------------------------------------------------------//
	class OBJnode
	{
	public:
		unsigned int index;
		bool averaged;
		OBJnode* next;

		OBJnode()
		{
			index = 0;
			next = NULL;
		}
		virtual ~OBJnode() {}
	};


	CAccessObj();
	virtual ~CAccessObj();

	COBJmodel *m_pModel;
	COBJgroup *m_pCurGroup;
	Vector3DF m_vMax, m_vMin;
	Vector3DF m_Res_TreGrid;

	void CalcBoundingBox();
	bool Equal(Vector3DF * u, Vector3DF * v, float epsilon);
	
	COBJgroup* FindGroup(char* name)
	{
		COBJgroup * group;
		
		assert(m_pModel);
		group = m_pModel->pGroups;
		while(group)
		{
			if (!strcmp(name, group->name))	break;
			group = group->next;
		}
		
		return group;
	}
	
	COBJgroup* AddGroup(char* name)
	{
		/*COBJgroup* group;
		group = FindGroup(name);
		if (!group)
		{
			group = new COBJgroup;
			sprintf(group->name, "%s", name);
			group->material = 0;
			group->nTriangles = 0;
			group->pTriangles = NULL;
			group->next = m_pModel->pGroups;
			m_pModel->pGroups = group;
			m_pModel->nGroups++;
		}*/

		COBJgroup* group;
		COBJgroup* tail;
		group = FindGroup(name);
		if (!group)
		{
			group = new COBJgroup;
			sprintf(group->name, "%s", name);
			group->material = -1;		//-1 means no material.
			group->nTriangles = 0;
			group->pTriangles = NULL;
			
			//Find the tail of group.
			tail = m_pModel->pGroups;
			if(tail == NULL)	//The model group is null.
				m_pModel->pGroups = group;
			else
			{
				while(tail->next != NULL)
					tail = tail->next;
				tail->next = group;
				group->next = NULL;
			}
			m_pModel->nGroups++;
		}
		
		return group;
	}
	
	unsigned int FindMaterial(char* name);
	char* DirName(char* path);
	void ReadMTL(char* name);
	void WriteMTL(char* modelpath, char* mtllibname);
	void FirstPass(FILE* file);
	void SecondPass(FILE* file);
	void Dimensions(float* dimensions);
	void Scale(float scale);
	void ReverseWinding();
	void FacetNormals();
	void VertexNormals(float angle);
	void LinearTexture();
	void SpheremapTexture();
	float Unitize();
	void WriteOBJ(char* filename, unsigned int mode);
	unsigned int OpenGLList();
	void WritePoly( char *filename );
	int GetGroupNumTri();
	unsigned int * GetGroupTriPIdx();
	void NextGroup();
	void FirstGroup();
	void Destroy();
	void Boundingbox(Vector3DF &vMax, Vector3DF &vMin);
	void CorrectMesh(float dis);

	void LoadOBJ(/*CString*/ char* filename/*, float minx, float miny, float minz, float maxx, float maxy, float maxz*/);
	void Draw();
	void DrawGrid();
	void setBoundingbox(Vector3DF max, Vector3DF min){ m_vMax = max; m_vMin = min; };
	void setCellSize(float size, Vector3DF delta){ CellSize = size;  gridDelta = delta; };

	void initTriGrid();
	void InsertGrid();
	int  ComputPos(Vector3DF pos, Vector3DF &Res);
	int  ComputPos(float *ipos, float *res);
	Vector3DF closest_Tri(int TriNum, Vector3DF p);
	float* closest_Tri(int TriNum, double* Pos);
	float  dist_Point2(Vector3DF a, Vector3DF b);
	bool   CohenSutherland(Vector3DF CellPos, Vector3DF p1, Vector3DF p2);
	void   fillGrid(Vector3DF p1, Vector3DF p2, Vector3DF pos1, Vector3DF pos2, int nTri);
	bool   IntersectSegmentTriangle(Vector3DF p, Vector3DF q, Vector3DF a, Vector3DF b, Vector3DF c,Vector3DF &Cp, float &t);
	bool   IntersectTest(Vector3DF prePos, Vector3DF pos, Vector3DF &n, Vector3DF &Cp,  float &t);
	bool   Collide_Tree(Vector3DF pos, Vector3DF &closeP, float &dis);
	Vector3DF computeDensity();
	Vector3DF computeRepAdh(Vector3DF pos, Vector3DF v,  float h, float scale, float space);
	Vector3DF computFric(Vector3DF pos, Vector3DF v, float h, float scale, float idensity, float m_LapKern, float  visc, float space);

public:
	int nCell;
	int* GridCell; // Vertices -> Grid;
	node** GridTri; //Grid -> Triangles; 一个三角形可能会和多个网格相交，所以不容易计算被占据的网格的个数， 哈希表，不好驾驭啊
	int m_GridAdj[30];
	float        Weight[8];
	Vector3DF    BarycCoord[8];
	float        WUnitMass;
	float        k , r0;
	int          m_GridAdjCnt;
	int          nadj;
	int*         findedTri;
	float        WoodDensity;
	float        CellSize;
	Vector3DF    gridDelta;   
};

#endif