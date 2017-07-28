//#include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "AccessObj.h"
#include <GL/glut.h>


#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CAccessObj::CAccessObj()
{
	m_pModel = NULL;

	GridTri = NULL;
	GridCell = NULL;

}

CAccessObj::~CAccessObj()
{
	Destroy();
}

//////////////////////////////////////////////////////////////////////
// Equal: compares two vectors and returns GL_TRUE if they are
// equal (within a certain threshold) or GL_FALSE if not. An epsilon
// that works fairly well is 0.000001.
//
// u - array of 3 GLfloats (float u[3])
// v - array of 3 GLfloats (float v[3]) 
//////////////////////////////////////////////////////////////////////
bool CAccessObj::Equal(Vector3DF * u, Vector3DF * v, float epsilon)
{
	if (objAbs(u->x - v->x) < epsilon &&
		objAbs(u->y - v->y) < epsilon &&
		objAbs(u->z - v->z) < epsilon) 
	{
		return GL_TRUE;
	}
	return GL_FALSE;
}


//////////////////////////////////////////////////////////////////////
// FindGroup: Find a group in the model
//////////////////////////////////////////////////////////////////////
/*COBJgroup * CAccessObj::FindGroup(char* name)
{
	COBJgroup * group;
	
	assert(m_pModel);
	
	group = m_pModel->pGroups;
	while(group) 
	{
		if (!strcmp(name, group->name))
			break;
		group = group->next;
	}
	
	return group;
}
*/
//////////////////////////////////////////////////////////////////////
// AddGroup: Add a group to the model
//////////////////////////////////////////////////////////////////////
/*COBJgroup * CAccessObj::AddGroup(char* name)
{
	COBJgroup* group;
	
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
	}
	
	return group;
}
*/
//////////////////////////////////////////////////////////////////////
// FindGroup: Find a material in the model
//////////////////////////////////////////////////////////////////////
unsigned int CAccessObj::FindMaterial(char* name)
{
	unsigned int i;
	bool bFound = false;
	
	// XXX doing a linear search on a string key'd list is pretty lame, but it works and is fast enough for now.
	for (i = 0; i < m_pModel->nMaterials; i++)
	{
		if (!strcmp(m_pModel->pMaterials[i].name, name))
		{
			bFound = true;
			break;
		}
	}
	
	// didn't find the name, so print a warning and return the default material (0).
	if (!bFound)
	{
		printf("FindMaterial():  can't find material \"%s\".\n", name);
		i = 0;
	}
	
	return i;
}


//////////////////////////////////////////////////////////////////////
// DirName: return the directory given a path
//
// path - filesystem path
//
// NOTE: the return value should be free'd.
//////////////////////////////////////////////////////////////////////
char * CAccessObj::DirName(char* path)
{
	char* dir;
	char* s;
	
	dir = strdup(path);
	
	s = strrchr(dir, '/');
	if (s)	s[1] = '\0';
	else	dir[0] = '\0';
	
	return dir;
}


//////////////////////////////////////////////////////////////////////
// ReadMTL: read a wavefront material library file
//
// model - properly initialized COBJmodel structure
// name  - name of the material library
//////////////////////////////////////////////////////////////////////
void CAccessObj::ReadMTL(char* name)
{
	FILE* file;
	char* dir;
	char* filename;
	char  buf[128];
	unsigned int nMaterials, i;
	
	dir = DirName(m_pModel->pathname);
	filename = new char [(strlen(dir) + strlen(name) + 1)];
	strcpy(filename, dir);
	strcat(filename, name);
	delete []dir;
	
	file = fopen(filename, "r");
	if (!file)
	{
		fprintf(stderr, "ReadMTL() failed: can't open material file \"%s\".\n", filename);
		exit(1);
	}
	delete [] filename;
	
	// count the number of materials in the file
	nMaterials = 1;
	while(fscanf(file, "%s", buf) != EOF)
	{
		switch(buf[0])
		{
		case '#':				/* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'n':				/* newmtl */
			fgets(buf, sizeof(buf), file);
			nMaterials++;
			sscanf(buf, "%s %s", buf, buf);
			break;
		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		}
	}
	
	rewind(file);
	
	m_pModel->pMaterials = new COBJmaterial [nMaterials];
	m_pModel->nMaterials = nMaterials;
	
	// set the default material
	for (i = 0; i < nMaterials; i++) {
		m_pModel->pMaterials[i].name[0] = '\0';
		m_pModel->pMaterials[i].shininess[0] = 65.0f;
		m_pModel->pMaterials[i].diffuse[0] = 0.8f;
		m_pModel->pMaterials[i].diffuse[1] = 0.8f;
		m_pModel->pMaterials[i].diffuse[2] = 0.8f;
		m_pModel->pMaterials[i].diffuse[3] = 1.0f;
		m_pModel->pMaterials[i].ambient[0] = 0.2f;
		m_pModel->pMaterials[i].ambient[1] = 0.2f;
		m_pModel->pMaterials[i].ambient[2] = 0.2f;
		m_pModel->pMaterials[i].ambient[3] = 1.0f;
		m_pModel->pMaterials[i].specular[0] = 0.0f;
		m_pModel->pMaterials[i].specular[1] = 0.0f;
		m_pModel->pMaterials[i].specular[2] = 0.0f;
		m_pModel->pMaterials[i].specular[3] = 1.0f;
		m_pModel->pMaterials[i].emissive[0] = 0.0f;
		m_pModel->pMaterials[i].emissive[1] = 0.0f;
		m_pModel->pMaterials[i].emissive[2] = 0.0f;
		m_pModel->pMaterials[i].emissive[3] = 0.0f;
	}
	sprintf(m_pModel->pMaterials[0].name, "default");
	
	// now, read in the data
	nMaterials = 0;
	while(fscanf(file, "%s", buf) != EOF)
	{
		switch(buf[0])
		{
		case '#':				/* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'n':				/* newmtl */
			fgets(buf, sizeof(buf), file);
			sscanf(buf, "%s %s", buf, buf);
			nMaterials++;
			sprintf(m_pModel->pMaterials[nMaterials].name, "%s", buf);
			break;
		case 'N':
			fscanf(file, "%f", *m_pModel->pMaterials[nMaterials].shininess);
			/* wavefront shininess is from [0, 1000], so scale for OpenGL */
			m_pModel->pMaterials[nMaterials].shininess[0] /= 1000.0;
			m_pModel->pMaterials[nMaterials].shininess[0] *= 128.0;
			break;

		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		}
	}
}

//////////////////////////////////////////////////////////////////////
// WriteMTL: write a wavefront material library file
//
// model      - properly initialized COBJmodel structure
// modelpath  - pathname of the model being written
// mtllibname - name of the material library to be written
//////////////////////////////////////////////////////////////////////
void CAccessObj::WriteMTL(char* modelpath, char* mtllibname)
{
	FILE* file;
	char* dir;
	char* filename;
	COBJmaterial* material;
	unsigned int i;
	
	dir = DirName(modelpath);
	filename = new char [(strlen(dir)+strlen(mtllibname))];
	strcpy(filename, dir);
	strcat(filename, mtllibname);
	delete []dir;
	
	/* open the file */
	file = fopen(filename, "w");
	if (!file)
	{
		fprintf(stderr, "WriteMTL() failed: can't open file \"%s\".\n", filename);
		exit(1);
	}
	delete [] filename;
	
	/* spit out a header */
	fprintf(file, "#  \n");
	fprintf(file, "#  Wavefront MTL generated by OBJ library\n");
	fprintf(file, "#  \n");
	fprintf(file, "#  OBJ library\n");
	fprintf(file, "#  Nate Robins\n");
	fprintf(file, "#  ndr@pobox.com\n");
	fprintf(file, "#  http://www.pobox.com/~ndr\n");
	fprintf(file, "#  \n\n");
	
	for (i = 0; i < m_pModel->nMaterials; i++)
	{
		material = &m_pModel->pMaterials[i];
		fprintf(file, "newmtl %s\n", material->name);
		fprintf(file, "Ka %f %f %f\n", 
			material->ambient[0], material->ambient[1], material->ambient[2]);
		fprintf(file, "Kd %f %f %f\n", 
			material->diffuse[0], material->diffuse[1], material->diffuse[2]);
		fprintf(file, "Ks %f %f %f\n", 
			material->specular[0],material->specular[1],material->specular[2]);
		fprintf(file, "Ns %f\n", material->shininess[0] / 128.0 * 1000.0);
		fprintf(file, "\n");
	}
}

//////////////////////////////////////////////////////////////////////
// FirstPass: first pass at a Wavefront OBJ file that gets all the
// statistics of the model (such as #vertices, #normals, etc)
//
// model - properly initialized COBJmodel structure
// file  - (fopen'd) file descriptor 
//////////////////////////////////////////////////////////////////////
void CAccessObj::FirstPass(FILE* file) 
{
	unsigned int    nVertices;		/* number of vertices in m_pModel */
	unsigned int    nNormals;		/* number of normals in m_pModel */
	unsigned int    nTexCoords;		/* number of texcoords in m_pModel */
	unsigned int    nTriangles;		/* number of triangles in m_pModel */
	COBJgroup* group;			    /* current group */
	unsigned  v, n, t;
	char      buf[128];
	
	/* make a default group */
	group = AddGroup("default");
	
	nVertices = nNormals = nTexCoords = nTriangles = 0;
	while(fscanf(file, "%s", buf) != EOF)
	{
		switch(buf[0])
		{
		case '#':				/* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'v':				/* v, vn, vt */
			switch(buf[1])
			{
			case '\0':			/* vertex */
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				nVertices++;
				break;
			case 'n':				/* normal */
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				nNormals++;
				break;
			case 't':				/* texcoord */
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				nTexCoords++;
				break;
			default:
				printf("FirstPass(): Unknown token \"%s\".\n", buf);
				exit(1);
				break;
			}
			break;

		case 'm':
			fgets(buf, sizeof(buf), file);
			sscanf(buf, "%s %s", buf, buf);
			sprintf(m_pModel->mtllibname, "%s", buf);
			ReadMTL(buf);
			break;
		case 'u':
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			m_pModel->m_nMaterial++;
			break;
		case 'g':				/* group */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			buf[strlen(buf)-1] = '\0';	/* nuke '\n' */
			group = AddGroup(buf);
			break;
		case 'f':				/* face */
			v = n = t = 0;
			fscanf(file, "%s", buf);
			/* can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d */
			if (strstr(buf, "//"))
			{
				/* v//n */
				sscanf(buf, "%d//%d", &v, &n);
				fscanf(file, "%d//%d", &v, &n);
				fscanf(file, "%d//%d", &v, &n);
				nTriangles++;
				group->nTriangles++;	//group中只记录面，不包含点.
				while(fscanf(file, "%d//%d", &v, &n) > 0)
				{
					nTriangles++;
					group->nTriangles++;
				}
			}
			else if (sscanf(buf, "%d/%d/%d", &v, &t, &n) == 3)
			{
				/* v/t/n */
				fscanf(file, "%d/%d/%d", &v, &t, &n);
				fscanf(file, "%d/%d/%d", &v, &t, &n);
				nTriangles++;
				group->nTriangles++;
				while(fscanf(file, "%d/%d/%d", &v, &t, &n) > 0)
				{
					nTriangles++;
					group->nTriangles++;
				}
			}
			else if (sscanf(buf, "%d/%d", &v, &t) == 2)
			{
				/* v/t */
				fscanf(file, "%d/%d", &v, &t);
				fscanf(file, "%d/%d", &v, &t);
				nTriangles++;
				group->nTriangles++;
				while(fscanf(file, "%d/%d", &v, &t) > 0)
				{
					nTriangles++;
					group->nTriangles++;
				}
			}
			else
			{
				/* v */
				fscanf(file, "%d", &v);
				fscanf(file, "%d", &v);
				nTriangles++;
				group->nTriangles++;
				while(fscanf(file, "%d", &v) > 0)
				{
					nTriangles++;
					group->nTriangles++;
				}
			}
			break;
			
		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		}
	}
	
	/* set the stats in the m_pModel structure */
	m_pModel->nVertices  = nVertices;
	m_pModel->nNormals   = nNormals;
	m_pModel->nTexCoords = nTexCoords;
	m_pModel->nTriangles = nTriangles;
	
	/* allocate memory for the triangles in each group */
	group = m_pModel->pGroups;
	while(group)
	{
		if (group->nTriangles)
			group->pTriangles = new unsigned int [group->nTriangles];
		group->nTriangles = 0;
		group = group->next;
	}
}

//////////////////////////////////////////////////////////////////////
// SecondPass: second pass at a Wavefront OBJ file that gets all
// the data.
//
// model - properly initialized COBJmodel structure
// file  - (fopen'd) file descriptor 
//////////////////////////////////////////////////////////////////////
void CAccessObj::SecondPass(FILE* file) 
{
	unsigned int	nVertices;		/* number of vertices in m_pModel */
	unsigned int	nNormals;		/* number of normals in m_pModel */
	unsigned int	nTexCoords;		/* number of texcoords in m_pModel */
	unsigned int	nTriangles;		/* number of triangles in m_pModel */
	Vector3DF *	vertices;		/* array of vertices  */
	Vector3DF *	normals;		/* array of normals */
	Vector3DF *	texcoords;		/* array of texture coordinates */
	COBJgroup *	group;			/* current group pointer */
	unsigned int	material;		/* current material */
	unsigned int	v, n, t;
	char		buf[128];
	
	/* set the pointer shortcuts */
	vertices     = m_pModel->vpVertices;
	normals      = m_pModel->vpNormals;
	texcoords    = m_pModel->vpTexCoords;
	group        = m_pModel->pGroups;
	
	/* on the second pass through the file, read all the data into the
	allocated arrays */
	nVertices = nNormals = nTexCoords = 1;
	nTriangles = 0;
	material = 0;

	while(fscanf(file, "%s", buf) != EOF)
	{
		switch(buf[0])
		{
		case '#':				/* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'v':				/* v, vn, vt */
			switch(buf[1])
			{
			case '\0':			/* vertex */
				fscanf(file, "%f %f %f", 
					&vertices[nVertices].x, 
					&vertices[nVertices].y, 
					&vertices[nVertices].z);
				nVertices++;
				break;
			case 'n':				/* normal */
				fscanf(file, "%f %f %f", 
					&normals[nNormals].x,
					&normals[nNormals].y, 
					&normals[nNormals].z);
				nNormals++;
				break;
			case 't':				/* texcoord */
				fscanf(file, "%f %f", 
					&texcoords[nTexCoords].x,
					&texcoords[nTexCoords].y);
				nTexCoords++;
				break;
			}
			break;
		case 'u':
			fgets(buf, sizeof(buf), file);
			sscanf(buf, "%s %s", buf, buf);
			group->material = material = FindMaterial(buf);
			break;
		case 'g':				/* group */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			buf[strlen(buf)-1] = '\0';	/* nuke '\n' */
			group = FindGroup(buf);
			group->material = material;
			break;
		case 'f':				/* face */
			v = n = t = 0;
			fscanf(file, "%s", buf);
			/* can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d */
			if (strstr(buf, "//"))
			{
				/* v//n */
				sscanf(buf, "%d//%d", &v, &n);
				T(nTriangles).VertIdx[0] = v;
				T(nTriangles).NormIdx[0] = n;
				fscanf(file, "%d//%d", &v, &n);
				T(nTriangles).VertIdx[1] = v;
				T(nTriangles).NormIdx[1] = n;
				fscanf(file, "%d//%d", &v, &n);
				T(nTriangles).VertIdx[2] = v;
				T(nTriangles).NormIdx[2] = n;
				group->pTriangles[group->nTriangles++] = nTriangles;
				nTriangles++;
				while(fscanf(file, "%d//%d", &v, &n) > 0)
				{
					T(nTriangles).VertIdx[0] = T(nTriangles-1).VertIdx[0];
					T(nTriangles).NormIdx[0] = T(nTriangles-1).NormIdx[0];
					T(nTriangles).VertIdx[1] = T(nTriangles-1).VertIdx[2];
					T(nTriangles).NormIdx[1] = T(nTriangles-1).NormIdx[2];
					T(nTriangles).VertIdx[2] = v;
					T(nTriangles).NormIdx[2] = n;
					group->pTriangles[group->nTriangles++] = nTriangles;
					nTriangles++;
				}
			}
			else if (sscanf(buf, "%d/%d/%d", &v, &t, &n) == 3)
			{
				/* v/t/n */
				T(nTriangles).VertIdx[0] = v;
				T(nTriangles).TexIdx[0] = t;
				T(nTriangles).NormIdx[0] = n;
				fscanf(file, "%d/%d/%d", &v, &t, &n);
				T(nTriangles).VertIdx[1] = v;
				T(nTriangles).TexIdx[1] = t;
				T(nTriangles).NormIdx[1] = n;
				fscanf(file, "%d/%d/%d", &v, &t, &n);
				T(nTriangles).VertIdx[2] = v;
				T(nTriangles).TexIdx[2] = t;
				T(nTriangles).NormIdx[2] = n;
				group->pTriangles[group->nTriangles++] = nTriangles;
				nTriangles++;
				while(fscanf(file, "%d/%d/%d", &v, &t, &n) > 0)
				{
					T(nTriangles).VertIdx[0] = T(nTriangles-1).VertIdx[0];
					T(nTriangles).TexIdx[0] = T(nTriangles-1).TexIdx[0];
					T(nTriangles).NormIdx[0] = T(nTriangles-1).NormIdx[0];
					T(nTriangles).VertIdx[1] = T(nTriangles-1).VertIdx[2];
					T(nTriangles).TexIdx[1] = T(nTriangles-1).TexIdx[2];
					T(nTriangles).NormIdx[1] = T(nTriangles-1).NormIdx[2];
					T(nTriangles).VertIdx[2] = v;
					T(nTriangles).TexIdx[2] = t;
					T(nTriangles).NormIdx[2] = n;
					group->pTriangles[group->nTriangles++] = nTriangles;
					nTriangles++;
				}
			}
			else if (sscanf(buf, "%d/%d", &v, &t) == 2)
			{
				/* v/t */
				T(nTriangles).VertIdx[0] = v;
				T(nTriangles).TexIdx[0] = t;
				fscanf(file, "%d/%d", &v, &t);
				T(nTriangles).VertIdx[1] = v;
				T(nTriangles).TexIdx[1] = t;
				fscanf(file, "%d/%d", &v, &t);
				T(nTriangles).VertIdx[2] = v;
				T(nTriangles).TexIdx[2] = t;
				group->pTriangles[group->nTriangles++] = nTriangles;
				nTriangles++;
				while(fscanf(file, "%d/%d", &v, &t) > 0)
				{
					T(nTriangles).VertIdx[0] = T(nTriangles-1).VertIdx[0];
					T(nTriangles).TexIdx[0] = T(nTriangles-1).TexIdx[0];
					T(nTriangles).VertIdx[1] = T(nTriangles-1).VertIdx[2];
					T(nTriangles).TexIdx[1] = T(nTriangles-1).TexIdx[2];
					T(nTriangles).VertIdx[2] = v;
					T(nTriangles).TexIdx[2] = t;
					group->pTriangles[group->nTriangles++] = nTriangles;
					nTriangles++;
				}
			}
			else
			{
				/* v */
				sscanf(buf, "%d", &v);
				T(nTriangles).VertIdx[0] = v;
				fscanf(file, "%d", &v);
				T(nTriangles).VertIdx[1] = v;
				fscanf(file, "%d", &v);
				T(nTriangles).VertIdx[2] = v;
				group->pTriangles[group->nTriangles++] = nTriangles;
				nTriangles++;
				while(fscanf(file, "%d", &v) > 0)
				{
					T(nTriangles).VertIdx[0] = T(nTriangles-1).VertIdx[0];
					T(nTriangles).VertIdx[1] = T(nTriangles-1).VertIdx[2];
					T(nTriangles).VertIdx[2] = v;
					group->pTriangles[group->nTriangles++] = nTriangles;
					nTriangles++;
				}
			}
			break;
			
		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		}
	}
}

/* public functions */

//////////////////////////////////////////////////////////////////////
// Unitize: "unitize" a model by translating it to the origin and
// scaling it to fit in a unit cube around the origin.  Returns the
// scalefactor used.
//
// model - properly initialized COBJmodel structure 
//////////////////////////////////////////////////////////////////////
float CAccessObj::Unitize()
{
 	unsigned int  i;
 	float maxx, minx, maxy, miny, maxz, minz;
 	float cx, cy, cz, w, h, d;
 	float scale;
 	
 	assert(m_pModel);
 	assert(m_pModel->vpVertices);
 	
 	/* get the max/mins */
 	maxx = minx = m_pModel->vpVertices[1].x;
 	maxy = miny = m_pModel->vpVertices[1].y;
	maxz = minz = m_pModel->vpVertices[1].z;
	for (i = 1; i <= m_pModel->nVertices; i++)
	/*maxx = minx = m_pModel->vpVertices[0].x;
 	maxy = miny = m_pModel->vpVertices[0].y;
	maxz = minz = m_pModel->vpVertices[0].z;
	for (i = 0; i < m_pModel->nVertices; i++)*/
	{
		if (maxx < m_pModel->vpVertices[i].x)
			maxx = m_pModel->vpVertices[i].x;
		if (minx > m_pModel->vpVertices[i].x)
			minx = m_pModel->vpVertices[i].x;
		
		if (maxy < m_pModel->vpVertices[i].y)
			maxy = m_pModel->vpVertices[i].y;
		if (miny > m_pModel->vpVertices[i].y)
			miny = m_pModel->vpVertices[i].y;
		
		if (maxz < m_pModel->vpVertices[i].z)
			maxz = m_pModel->vpVertices[i].z;
		if (minz > m_pModel->vpVertices[i].z)
			minz = m_pModel->vpVertices[i].z;
	}
	
	/* calculate m_pModel width, height, and depth */
	w = maxx - minx;
	h = maxy-miny;
	d = maxz-minz;
	
	/* calculate center of the m_pModel */
	cx = (maxx + minx) / 2.0f;
	cy = -5/*(maxy + miny) / 2.0f*/;
	cz = (maxz + minz) / 2.0f;
	
	/* calculate unitizing scale factor */
	scale = 2;//1.0f / objMax(objMax(w, h), d);
	
	/* translate around center then scale */
	for (i = 1; i <= m_pModel->nVertices; i++)
	//for (i = 0; i < m_pModel->nVertices; i++)
	{
		m_pModel->vpVertices[i].x -= cx;
		m_pModel->vpVertices[i].y -= cy;
		m_pModel->vpVertices[i].z -= cz;
		m_pModel->vpVertices[i].x *= scale;
		m_pModel->vpVertices[i].y *= scale;
		m_pModel->vpVertices[i].z *= scale;
	}
	
	return scale;
}


//////////////////////////////////////////////////////////////////////
// Dimensions: Calculates the dimensions (width, height, depth) of
// a model.
//
// model      - initialized COBJmodel structure
// dimensions - array of 3 GLfloats (float dimensions[3])
//////////////////////////////////////////////////////////////////////
void CAccessObj::Dimensions(float* dimensions)
{
	unsigned int i;
	Vector3DF vMax, vMin;
	
	assert(m_pModel);
	assert(m_pModel->vpVertices);
	assert(dimensions);
	
	/* get the max/mins */
	vMax = vMin = m_pModel->vpVertices[0];
	for (i = 1; i <= m_pModel->nVertices; i++)
	//for (i = 0; i < m_pModel->nVertices; i++)
	{
		if (vMax.x < m_pModel->vpVertices[i].x)
			vMax.x = m_pModel->vpVertices[i].x;
		if (vMin.x > m_pModel->vpVertices[i].x)
			vMin.x = m_pModel->vpVertices[i].x;
		
		if (vMax.y < m_pModel->vpVertices[i].y)
			vMax.y = m_pModel->vpVertices[i].y;
		if (vMin.y > m_pModel->vpVertices[i].y)
			vMin.y = m_pModel->vpVertices[i].y;
		
		if (vMax.z < m_pModel->vpVertices[i].z)
			vMax.z = m_pModel->vpVertices[i].z;
		if (vMin.z > m_pModel->vpVertices[i].z)
			vMin.z = m_pModel->vpVertices[i].z;
	}
	
	/* calculate m_pModel width, height, and depth */
	dimensions[0] = vMax.x-vMin.x;
	dimensions[1] = vMax.y-vMin.y;
	dimensions[2] = vMax.z-vMin.z;
}

//////////////////////////////////////////////////////////////////////
// Scale: Scales a model by a given amount.
// 
// model - properly initialized COBJmodel structure
// scale - scalefactor (0.5 = half as large, 2.0 = twice as large)
//////////////////////////////////////////////////////////////////////
void CAccessObj::Scale(float scale)
{
	unsigned int i;
	
	for (i = 1; i <= m_pModel->nVertices; i++)
	//for (i = 0; i < m_pModel->nVertices; i++)
		m_pModel->vpVertices[i] = m_pModel->vpVertices[i] * scale;
}

//////////////////////////////////////////////////////////////////////
// ReverseWinding: Reverse the polygon winding for all polygons in
// this model.  Default winding is counter-clockwise.  Also changes
// the direction of the normals.
// 
// model - properly initialized COBJmodel structure 
//////////////////////////////////////////////////////////////////////
void CAccessObj::ReverseWinding()
{
	unsigned int i, swap;
	
	assert(m_pModel);
	
	for (i = 0; i < m_pModel->nTriangles; i++)
	{
		swap = T(i).VertIdx[0];
		T(i).VertIdx[0] = T(i).VertIdx[2];
		T(i).VertIdx[2] = swap;
		
		if (m_pModel->nNormals)
		{
			swap = T(i).NormIdx[0];
			T(i).NormIdx[0] = T(i).NormIdx[2];
			T(i).NormIdx[2] = swap;
		}
		
		if (m_pModel->nTexCoords)
		{
			swap = T(i).TexIdx[0];
			T(i).TexIdx[0] = T(i).TexIdx[2];
			T(i).TexIdx[2] = swap;
		}
	}
	
	/* reverse facet normals */
	for (i = 1; i <= m_pModel->nFacetnorms; i++)
		m_pModel->vpFacetNorms[i] = m_pModel->vpFacetNorms[i]*(-1);
	
	/* reverse vertex normals */
	for (i = 1; i <= m_pModel->nNormals; i++)
		//for (i = 0; i < m_pModel->nNormals; i++)
		m_pModel->vpNormals[i] = m_pModel->vpNormals[i]*(-1);
}

//////////////////////////////////////////////////////////////////////
// FacetNormals: Generates facet normals for a model (by taking the
// cross product of the two vectors derived from the sides of each
// triangle).  Assumes a counter-clockwise winding.
//
// model - initialized COBJmodel structure
//////////////////////////////////////////////////////////////////////
void CAccessObj::FacetNormals()
{
	unsigned int  i;
	Vector3DF u, v;
	
	assert(m_pModel);
	assert(m_pModel->vpVertices);
	
	/* clobber any old facetnormals */
	if (m_pModel->vpFacetNorms)
		delete[] m_pModel->vpFacetNorms;
	
	/* allocate memory for the new facet normals */
	m_pModel->nFacetnorms = m_pModel->nTriangles;
	m_pModel->vpFacetNorms = new Vector3DF [m_pModel->nFacetnorms + 1];
	
	for (i = 0; i < m_pModel->nTriangles; i++)
	{
		m_pModel->pTriangles[i].findex = i+1;
		
		u = m_pModel->vpVertices[T(i).VertIdx[1]] - m_pModel->vpVertices[T(i).VertIdx[0]];
		v = m_pModel->vpVertices[T(i).VertIdx[2]] - m_pModel->vpVertices[T(i).VertIdx[0]];

		m_pModel->vpFacetNorms[i+1] = u.Cross(v);
		m_pModel->vpFacetNorms[i+1].Normalize();
	}
}

//////////////////////////////////////////////////////////////////////
// VertexNormals: Generates smooth vertex normals for a model.
// First builds a list of all the triangles each vertex is in.  Then
// loops through each vertex in the the list averaging all the facet
// normals of the triangles each vertex is in.  Finally, sets the
// normal index in the triangle for the vertex to the generated smooth
// normal. 
// If the dot product of a facet normal and the facet normal
// associated with the 【first triangle in the list of triangles the
// current vertex is in】 is greater than the cosine of the angle
// parameter to the function, that facet normal is not added into the
// average normal calculation and the corresponding vertex is given
// the facet normal.  This tends to preserve hard edges.  The angle to
// use depends on the model, but 90 degrees is usually a good start.
//
// model - initialized COBJmodel structure
// angle - maximum angle (in degrees) to smooth across
//////////////////////////////////////////////////////////////////////
void CAccessObj::VertexNormals(float angle)
{
	OBJnode*  node;
	OBJnode*  tail;
	OBJnode** members;
	Vector3DF *  normals;
	unsigned int    nNormals;
	Vector3DF   average;
	float   dot, cos_angle;
	unsigned int    i, avg;
	
	assert(m_pModel);
	assert(m_pModel->vpFacetNorms);
	
	/* calculate the cosine of the angle (in degrees) */
	cos_angle = (float)cos(angle * 3.14159265f / 180.0f);
	
	/* nuke any previous normals */
	if (m_pModel->vpNormals)
		delete []m_pModel->vpNormals;
	
	/* allocate space for new normals */
	m_pModel->nNormals = m_pModel->nTriangles * 3; /* 3 normals per triangle */
	m_pModel->vpNormals = new Vector3DF [m_pModel->nNormals+1];
	//m_pModel->vpNormals = new Vector3DF [m_pModel->nNormals];
	
	/* allocate a structure that will hold a linked list of triangle
	indices for each vertex */
	members = new OBJnode * [m_pModel->nVertices + 1];
	for (i = 1; i <= m_pModel->nVertices; i++)
	//members = new OBJnode * [m_pModel->nVertices];
	//for (i = 0; i < m_pModel->nVertices; i++)
		members[i] = NULL;
	
	/* for every triangle, create a node for each vertex in it */
	for (i = 0; i < m_pModel->nTriangles; i++)//对每个三角形，都创建三个顶点的节点并插入链表
	{
		node = new OBJnode;//声明一个节点指向members中的顶点前面，即一个数组里都是链表头
		node->index = i;
		node->next  = members[T(i).VertIdx[0]];
		members[T(i).VertIdx[0]] = node;
		
		node = new OBJnode;
		node->index = i;
		node->next  = members[T(i).VertIdx[1]];
		members[T(i).VertIdx[1]] = node;
		
		node = new OBJnode;
		node->index = i;
		node->next  = members[T(i).VertIdx[2]];
		members[T(i).VertIdx[2]] = node;
	}
	
	/* calculate the average normal for each vertex */
	nNormals = 1;
	for (i = 1; i <= m_pModel->nVertices; i++)
	//nNormals = 0;
	//for (i = 0; i < m_pModel->nVertices; i++)
	{
		// calculate an average normal for this vertex by averaging the
		// facet normal of every triangle this vertex is in
		node = members[i];
		if (!node)
			fprintf(stderr, "VertexNormals(): vertex w/o a triangle\n");
		average = Vector3DF(0, 0, 0);
		avg = 0;
		while (node)//第一遍循环求平均
		{
			/* only average if the dot product of the angle between the two
			facet normals is greater than the cosine of the threshold
			angle -- or, said another way, the angle between the two
			facet normals is less than (or equal to) the threshold angle */
			dot = m_pModel->vpFacetNorms[T(node->index).findex].Dot( m_pModel->vpFacetNorms[T(members[i]->index).findex]);
			if (dot > cos_angle)
			{
				node->averaged = GL_TRUE;
				average = average + m_pModel->vpFacetNorms[T(node->index).findex];
				avg = 1;			/* we averaged at least one normal! */
			}
			else
			{
				node->averaged = GL_FALSE;
			}
			node = node->next;
		}
		
		if (avg)
		{
			/* normalize the averaged normal */
			average.Normalize();
			
			/* add the normal to the vertex normals list */
			m_pModel->vpNormals[nNormals] = average;
			avg = nNormals;
			nNormals++;
		}
		
		/* set the normal of this vertex in each triangle it is in */
		node = members[i];
		while (node)
		{
			if (node->averaged)//更新三角面片中的顶点的向量
			{
				/* if this node was averaged, use the average normal */
				if (T(node->index).VertIdx[0] == i)
					T(node->index).NormIdx[0] = avg;
				else if (T(node->index).VertIdx[1] == i)
					T(node->index).NormIdx[1] = avg;
				else if (T(node->index).VertIdx[2] == i)
					T(node->index).NormIdx[2] = avg;
			}
			else
			{
				/* if this node wasn't averaged, use the facet normal */
				m_pModel->vpNormals[nNormals] = m_pModel->vpFacetNorms[T(node->index).findex];
				if (T(node->index).VertIdx[0] == i)
					T(node->index).NormIdx[0] = nNormals;
				else if (T(node->index).VertIdx[1] == i)
					T(node->index).NormIdx[1] = nNormals;
				else if (T(node->index).VertIdx[2] == i)
					T(node->index).NormIdx[2] = nNormals;
				nNormals++;
			}
			node = node->next;
		}
	}
	
	m_pModel->nNormals = nNormals - 1;
	//m_pModel->nNormals = nNormals;
	
	/* free the member information */
	for (i = 1; i <= m_pModel->nVertices; i++)
	//for (i = 0; i < m_pModel->nVertices; i++)
	{
		node = members[i];
		while (node)
		{
			tail = node;
			node = node->next;
			delete tail;
		}
	}
	delete []members;
	
	/* pack the normals array (we previously allocated the maximum
	number of normals that could possibly be created (nTriangles *
	3), so get rid of some of them (usually alot unless none of the
	facet normals were averaged)) */
	normals = m_pModel->vpNormals;
	m_pModel->vpNormals = new Vector3DF [m_pModel->nNormals+1];
	for (i = 1; i <= m_pModel->nNormals; i++) {
	/*m_pModel->vpNormals = new Vector3DF [m_pModel->nNormals];
	for (i = 0; i < m_pModel->nNormals; i++) {*/
		m_pModel->vpNormals[i] = normals[i];
	}
	delete []normals;
}


//////////////////////////////////////////////////////////////////////
// LinearTexture: Generates texture coordinates according to a
// linear projection of the texture map.  It generates these by
// linearly mapping the vertices onto a square.
//
// model - pointer to initialized COBJmodel structure
//////////////////////////////////////////////////////////////////////
void CAccessObj::LinearTexture()
{
	COBJgroup *group;
	float dimensions[3];
	float x, y, scalefactor;
	unsigned int i;
	
	assert(m_pModel);
	
//	if (m_pModel->vpTexCoords)
//		free(m_pModel->vpTexCoords);
//	m_pModel->nTexCoords = m_pModel->nVertices;
//	m_pModel->vpTexCoords = new Vector3DF [m_pModel->nTexCoords+1];
	
	//if (m_pModel->vpTexCoords)
	//	free(m_pModel->vpTexCoords);
	//m_pModel->nTexCoords = m_pModel->nVertices;
	//m_pModel->vpTexCoords = new Vector3DF [m_pModel->nTexCoords];

	Dimensions(dimensions);
	scalefactor = 2.0f / objAbs(objMax(objMax(dimensions[0], dimensions[1]), dimensions[2]));
	
	/* do the calculations */
	for(i = 1; i <= m_pModel->nVertices; i++)
	//for(i = 0; i < m_pModel->nVertices; i++)
	{
		x = m_pModel->vpVertices[i].x * scalefactor;
		y = m_pModel->vpVertices[i].z * scalefactor;
//		m_pModel->vpTexCoords[i] = Vector3DF((x + 1.0) / 2.0, (y + 1.0) / 2.0, 0.0);
		//m_pModel->vpTexCoords[i] = Vector3DF((x + 1.0) / 2.0, (y + 1.0) / 2.0, 0.0);
	}
	
	/* go through and put texture coordinate indices in all the triangles */
	group = m_pModel->pGroups;
	while(group)
	{
		for(i = 0; i < group->nTriangles; i++)
		{
			T(group->pTriangles[i]).TexIdx[0] = T(group->pTriangles[i]).VertIdx[0];
			T(group->pTriangles[i]).TexIdx[1] = T(group->pTriangles[i]).VertIdx[1];
			T(group->pTriangles[i]).TexIdx[2] = T(group->pTriangles[i]).VertIdx[2];
		}    
		group = group->next;
	}	
}

//////////////////////////////////////////////////////////////////////
// SpheremapTexture: Generates texture coordinates according to a
// spherical projection of the texture map.  Sometimes referred to as
// spheremap, or reflection map texture coordinates.  It generates
// these by using the normal to calculate where that vertex would map
// onto a sphere.  Since it is impossible to map something flat
// perfectly onto something spherical, there is distortion at the
// poles.  This particular implementation causes the poles along the X
// axis to be distorted.
//
// model - pointer to initialized COBJmodel structure
//////////////////////////////////////////////////////////////////////
void CAccessObj::SpheremapTexture()
{
	COBJgroup* group;
	float theta, phi, rho, x, y, z, r;
	unsigned int i;
	
	assert(m_pModel);
	assert(m_pModel->vpNormals);
	
//	if (m_pModel->vpTexCoords)
//		free(m_pModel->vpTexCoords);
//	m_pModel->nTexCoords = m_pModel->nNormals;
//	m_pModel->vpTexCoords = new Vector3DF [m_pModel->nTexCoords+1];

	if (m_pModel->vpTexCoords)
		free(m_pModel->vpTexCoords);
	m_pModel->nTexCoords = m_pModel->nNormals;
	m_pModel->vpTexCoords = new Vector3DF [m_pModel->nTexCoords+1];
	//m_pModel->vpTexCoords = new Vector3DF [m_pModel->nTexCoords];
	
	//for (i = 0; i < m_pModel->nNormals; i++)
	for (i = 1; i <= m_pModel->nNormals; i++)
	{
		z = m_pModel->vpNormals[i].x;	/* re-arrange for pole distortion */
		y = m_pModel->vpNormals[i].y;
		x = m_pModel->vpNormals[i].z;
		r = (float)sqrt((x * x) + (y * y));
		rho = (float)sqrt((r * r) + (z * z));
		
		if(r == 0.0)
		{
			theta = 0.0;
			phi = 0.0;
		}
		else
		{
			if(z == 0.0)
				phi = 3.14159265f / 2.0f;
			else	phi = (float)acos(z / rho);
			
			if(y == 0.0)
				theta = 3.141592365f / 2.0f;
			else	theta = (float)asin(y / r) + (3.14159f / 2.0f);
		}
		
//		m_pModel->vpTexCoords[i] = Vector3DF(theta / 3.14159f, phi / 3.14159265f, 0.0f);
	}
	
	/* go through and put texcoord indices in all the triangles */
	group = m_pModel->pGroups;
	while(group)
	{
		for (i = 0; i < group->nTriangles; i++)
		{
			T(group->pTriangles[i]).TexIdx[0] = T(group->pTriangles[i]).NormIdx[0];
			T(group->pTriangles[i]).TexIdx[1] = T(group->pTriangles[i]).NormIdx[1];
			T(group->pTriangles[i]).TexIdx[2] = T(group->pTriangles[i]).NormIdx[2];
		}
		group = group->next;
	}
}

//////////////////////////////////////////////////////////////////////
// objDelete: Deletes a COBJmodel structure.
//
// model - initialized COBJmodel structure
//////////////////////////////////////////////////////////////////////
void CAccessObj::Destroy()
{
	delete m_pModel;
	m_pModel = NULL;
}

//////////////////////////////////////////////////////////////////////
// objReadOBJ: Reads a model description from a Wavefront .OBJ file.
// Returns a pointer to the created object which should be free'd with
// objDelete().
//
// filename - name of the file containing the Wavefront .OBJ format data.  
//////////////////////////////////////////////////////////////////////
void CAccessObj::LoadOBJ(/*CString*/ char*  filename/*, float minx, float miny, float minz, float maxx, float maxy, float maxz*/)
{
	FILE*     file;
	
	// open the file
	file = fopen(filename, "r");
	if (!file)
	{
		fprintf(stderr, "objReadOBJ() failed: can't open data file \"%s\".\n",
			filename);
		exit(1);
	}
	
	if (m_pModel != NULL) Destroy();

	// allocate a new model
	m_pModel = new COBJmodel;

	sprintf(m_pModel->pathname, "%s", filename);
	m_pModel->mtllibname[0] = '\0';

	m_pModel->nVertices	= 0;
	m_pModel->vpVertices	= NULL;
	m_pModel->nNormals	= 0;
	m_pModel->vpNormals	= NULL;
	m_pModel->nTexCoords	= 0;
	m_pModel->vpTexCoords	= NULL;
	m_pModel->nFacetnorms	= 0;
	m_pModel->vpFacetNorms	= NULL;
	m_pModel->nTriangles	= 0;
	m_pModel->pTriangles	= NULL;
	m_pModel->nMaterials	= 0;
	m_pModel->pMaterials	= NULL;
	m_pModel->nGroups	= 0;
	m_pModel->pGroups	= NULL;
	m_pModel->position	= Vector3DF (0, 0, 0);
	
	// make a first pass through the file to get a count of the number
	// of vertices, normals, texcoords * triangles
	FirstPass(file);
	
	/* allocate memory */
	m_pModel->vpVertices = new Vector3DF [m_pModel->nVertices + 1];
	//m_pModel->vpVertices = new Vector3DF [m_pModel->nVertices];
	m_pModel->pTriangles = new COBJtriangle [m_pModel->nTriangles];
	if (m_pModel->nNormals)
	{
		m_pModel->vpNormals = new Vector3DF [m_pModel->nNormals + 1];
		//m_pModel->vpNormals = new Vector3DF [m_pModel->nNormals];
	}
	if (m_pModel->nTexCoords)
	{
		m_pModel->vpTexCoords = new Vector3DF [m_pModel->nTexCoords + 1];
		//m_pModel->vpTexCoords = new Vector3DF [m_pModel->nTexCoords];
	}
	
	/* rewind to beginning of file and read in the data this pass */
	rewind(file);
	SecondPass(file);
	
	/* close the file */
	fclose(file);

	Unitize();

////////////My add//////////////////////////////////////////////////
	FacetNormals();
	VertexNormals(90);
////////////////////////////////////////////////////////////////////

	//CalcBoundingBox();
	/*if(m_vMax.x < maxx) m_vMax.x = maxx;
	if(m_vMax.y < maxy) m_vMax.y = maxy;
	if(m_vMax.z < maxz) m_vMax.z = maxz;
	if(m_vMin.x > minx) m_vMin.x = minx;
	if(m_vMin.y > miny) m_vMin.y = miny;
	if(m_vMin.z > minz) m_vMin.z = minz;*/

	initTriGrid();
	InsertGrid();

}

//////////////////////////////////////////////////////////////////////
// WriteOBJ: Writes a model description in Wavefront .OBJ format to
// a file.
//
// model    - initialized COBJmodel structure
// filename - name of the file to write the Wavefront .OBJ format data to
// mode     - a bitwise or of values describing what is written to the file
//            OBJ_NONE     -  render with only vertices
//            OBJ_FLAT     -  render with facet normals
//            OBJ_SMOOTH   -  render with vertex normals
//            OBJ_TEXTURE  -  render with texture coords
//            OBJ_COLOR    -  render with colors (color material)
//            OBJ_MATERIAL -  render with materials
//            OBJ_COLOR and OBJ_MATERIAL should not both be specified.  
//            OBJ_FLAT and OBJ_SMOOTH should not both be specified.  
//////////////////////////////////////////////////////////////////////
void CAccessObj::WriteOBJ(char* filename, unsigned int mode)
{
	unsigned int    i;
	FILE*     file;
	COBJgroup* group;
	
	assert(m_pModel);
	
	/* do a bit of warning */
	if (mode * OBJ_FLAT && !m_pModel->vpFacetNorms)
	{
		printf("WriteOBJ() warning: flat normal output requested "
			"with no facet normals defined.\n");
		mode *= ~OBJ_FLAT;
	}
	if (mode * OBJ_SMOOTH && !m_pModel->vpNormals)
	{
		printf("WriteOBJ() warning: smooth normal output requested "
			"with no normals defined.\n");
		mode *= ~OBJ_SMOOTH;
	}
//	if (mode * OBJ_TEXTURE && !m_pModel->vpTexCoords)
//	{
//		printf("WriteOBJ() warning: texture coordinate output requested "
//			"with no texture coordinates defined.\n");
//		mode *= ~OBJ_TEXTURE;
//	}
	if (mode * OBJ_FLAT && mode * OBJ_SMOOTH)
	{
		printf("WriteOBJ() warning: flat normal output requested "
			"and smooth normal output requested (using smooth).\n");
		mode *= ~OBJ_FLAT;
	}
	if (mode * OBJ_COLOR && !m_pModel->pMaterials)
	{
		printf("WriteOBJ() warning: color output requested "
			"with no colors (materials) defined.\n");
		mode *= ~OBJ_COLOR;
	}
	if (mode * OBJ_MATERIAL && !m_pModel->pMaterials)
	{
		printf("WriteOBJ() warning: material output requested "
			"with no materials defined.\n");
		mode *= ~OBJ_MATERIAL;
	}
	if (mode * OBJ_COLOR && mode * OBJ_MATERIAL)
	{
		printf("WriteOBJ() warning: color and material output requested "
			"outputting only materials.\n");
		mode *= ~OBJ_COLOR;
	}
	
	
	/* open the file */
	file = fopen(filename, "w");
	if (!file)
	{
		fprintf(stderr, "WriteOBJ() failed: can't open file \"%s\" to write.\n",
			filename);
		exit(1);
	}
	
	/* spit out a header */
	fprintf(file, "#  \n");
	fprintf(file, "#  Wavefront OBJ\n");
	fprintf(file, "#  \n");
	
	if (mode * OBJ_MATERIAL && m_pModel->mtllibname)
	{
		fprintf(file, "\nmtllib %s\n\n", m_pModel->mtllibname);
		WriteMTL(filename, m_pModel->mtllibname);
	}
	
	/* spit out the vertices */
	fprintf(file, "\n");
	fprintf(file, "# %d vertices\n", m_pModel->nVertices);
	for (i = 1; i <= m_pModel->nVertices; i++)
	//for (i = 0; i < m_pModel->nVertices; i++)
	{
		fprintf(file, "v %f %f %f\n", 
			m_pModel->vpVertices[i].x,
			m_pModel->vpVertices[i].y,
			m_pModel->vpVertices[i].z);
	}
	
	/* spit out the smooth/flat normals */
	if (mode * OBJ_SMOOTH)
	{
		fprintf(file, "\n");
		fprintf(file, "# %d normals\n", m_pModel->nNormals);
		for (i = 1; i <= m_pModel->nNormals; i++)
		//for (i = 0; i < m_pModel->nNormals; i++)
		{
			fprintf(file, "vn %f %f %f\n", 
				m_pModel->vpNormals[i].x,
				m_pModel->vpNormals[i].y,
				m_pModel->vpNormals[i].z);
		}
	}
	else if (mode * OBJ_FLAT)
	{
		fprintf(file, "\n");
		fprintf(file, "# %d normals\n", m_pModel->nFacetnorms);
		for (i = 1; i <= m_pModel->nNormals; i++)
		//for (i = 0; i < m_pModel->nNormals; i++)
		{
			fprintf(file, "vn %f %f %f\n", 
				m_pModel->vpFacetNorms[i].x,
				m_pModel->vpFacetNorms[i].y,
				m_pModel->vpFacetNorms[i].z);
		}
	}
	
	/* spit out the texture coordinates */

	//if (mode * OBJ_TEXTURE)
	//{
	//	fprintf(file, "\n");
	//	fprintf(file, "# %d texcoords\n", m_pModel->vpTexCoords);
	//	//for (i = 1; i <= m_pModel->nTexCoords; i++)
	//	for (i = 0; i < m_pModel->nTexCoords; i++)
	//	{
	//		fprintf(file, "vt %f %f\n", 
	//			m_pModel->vpTexCoords[i].x,
	//			m_pModel->vpTexCoords[i].y);
	//	}
	//}

	
	fprintf(file, "\n");
	fprintf(file, "# %d groups\n", m_pModel->nGroups);
	fprintf(file, "# %d faces (triangles)\n", m_pModel->nTriangles);
	fprintf(file, "\n");
	
	group = m_pModel->pGroups;
	while(group)
	{
		fprintf(file, "g %s\n", group->name);
		if (mode * OBJ_MATERIAL)
			fprintf(file, "usemtl %s\n", m_pModel->pMaterials[group->material].name);
		for (i = 0; i < group->nTriangles; i++)
		{
			if (mode * OBJ_SMOOTH && mode * OBJ_TEXTURE)
			{
				fprintf(file, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
					T(group->pTriangles[i]).VertIdx[0], 
					T(group->pTriangles[i]).NormIdx[0], 
					T(group->pTriangles[i]).TexIdx[0],
					T(group->pTriangles[i]).VertIdx[1],
					T(group->pTriangles[i]).NormIdx[1],
					T(group->pTriangles[i]).TexIdx[1],
					T(group->pTriangles[i]).VertIdx[2],
					T(group->pTriangles[i]).NormIdx[2],
					T(group->pTriangles[i]).TexIdx[2]);
			}
			else if (mode * OBJ_FLAT && mode * OBJ_TEXTURE)
			{
				fprintf(file, "f %d/%d %d/%d %d/%d\n",
					T(group->pTriangles[i]).VertIdx[0],
					T(group->pTriangles[i]).findex,
					T(group->pTriangles[i]).VertIdx[1],
					T(group->pTriangles[i]).findex,
					T(group->pTriangles[i]).VertIdx[2],
					T(group->pTriangles[i]).findex);
			}
			else if (mode * OBJ_TEXTURE)
			{
				fprintf(file, "f %d/%d %d/%d %d/%d\n",
					T(group->pTriangles[i]).VertIdx[0],
					T(group->pTriangles[i]).TexIdx[0],
					T(group->pTriangles[i]).VertIdx[1],
					T(group->pTriangles[i]).TexIdx[1],
					T(group->pTriangles[i]).VertIdx[2],
					T(group->pTriangles[i]).TexIdx[2]);
			}
			else if (mode * OBJ_SMOOTH)
			{
				fprintf(file, "f %d//%d %d//%d %d//%d\n",
					T(group->pTriangles[i]).VertIdx[0],
					T(group->pTriangles[i]).NormIdx[0],
					T(group->pTriangles[i]).VertIdx[1],
					T(group->pTriangles[i]).NormIdx[1],
					T(group->pTriangles[i]).VertIdx[2], 
					T(group->pTriangles[i]).NormIdx[2]);
			}
			else if (mode * OBJ_FLAT)
			{
				fprintf(file, "f %d//%d %d//%d %d//%d\n",
					T(group->pTriangles[i]).VertIdx[0], 
					T(group->pTriangles[i]).findex,
					T(group->pTriangles[i]).VertIdx[1],
					T(group->pTriangles[i]).findex,
					T(group->pTriangles[i]).VertIdx[2],
					T(group->pTriangles[i]).findex);
			}
			else
			{
				fprintf(file, "f %d %d %d\n",
					T(group->pTriangles[i]).VertIdx[0],
					T(group->pTriangles[i]).VertIdx[1],
					T(group->pTriangles[i]).VertIdx[2]);
			}
		}
		fprintf(file, "\n");
		group = group->next;
	}
	
	fclose(file);
}

//////////////////////////////////////////////////////////////////////
void CAccessObj::Draw()
{
	static unsigned int i;
	static COBJgroup* group;
	static COBJtriangle* triangle;
	static COBJmaterial* material;

	if (m_pModel == NULL) return;
	//DrawGrid();
	group = m_pModel->pGroups;
	glDisable(GL_LIGHTING);
	glColor3f(0.0f, 1.0f, 0.3f);
	while (group)
	{
		for (i = 0; i < group->nTriangles; i++)
		{
			triangle = &T(group->pTriangles[i]);
			
		glBegin(/*GL_LINE_LOOP*/GL_TRIANGLES);
			glVertex3f(m_pModel->vpVertices[triangle->VertIdx[0]].x,
				m_pModel->vpVertices[triangle->VertIdx[0]].y,
				m_pModel->vpVertices[triangle->VertIdx[0]].z);

			glVertex3f(m_pModel->vpVertices[triangle->VertIdx[1]].x,
				m_pModel->vpVertices[triangle->VertIdx[1]].y,
				m_pModel->vpVertices[triangle->VertIdx[1]].z);

			glVertex3f(m_pModel->vpVertices[triangle->VertIdx[2]].x,
				m_pModel->vpVertices[triangle->VertIdx[2]].y,
				m_pModel->vpVertices[triangle->VertIdx[2]].z);
		glEnd();
		}
		
		group = group->next;
	}

	
}
void CAccessObj::DrawSmooth(){
	static unsigned int i;
	static COBJgroup* group;
	static COBJtriangle* triangle;
	static COBJmaterial* material;

	if (m_pModel == NULL) return;
	//DrawGrid();
	group = m_pModel->pGroups;
	glPushMatrix();
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	float diff[4] = {0.0, 0.8, 0, 1};
	glMaterialfv(GL_FRONT,GL_DIFFUSE,diff); 
	glColor3f(0.0f, 1.0f, 0.3f);
	while (group)
	{
		for (i = 0; i < group->nTriangles; i++)
		{
			triangle = &T(group->pTriangles[i]);

			glBegin(/*GL_LINE_LOOP*/GL_TRIANGLES);
			//glNormal3f(m_pModel->vpFacetNorms[triangle->findex].x, m_pModel->vpFacetNorms[triangle->findex].y, m_pModel->vpFacetNorms[triangle->findex].z);
			
			glNormal3f(m_pModel->vpNormals[triangle->NormIdx[0]].x, 
				m_pModel->vpNormals[triangle->NormIdx[0]].y,
				m_pModel->vpNormals[triangle->NormIdx[0]].z);
			glVertex3f(m_pModel->vpVertices[triangle->VertIdx[0]].x,
				m_pModel->vpVertices[triangle->VertIdx[0]].y,
				m_pModel->vpVertices[triangle->VertIdx[0]].z);

			glNormal3f(m_pModel->vpNormals[triangle->NormIdx[1]].x, 
				m_pModel->vpNormals[triangle->NormIdx[1]].y,
				m_pModel->vpNormals[triangle->NormIdx[1]].z);
			glVertex3f(m_pModel->vpVertices[triangle->VertIdx[1]].x,
				m_pModel->vpVertices[triangle->VertIdx[1]].y,
				m_pModel->vpVertices[triangle->VertIdx[1]].z);

			glNormal3f(m_pModel->vpNormals[triangle->NormIdx[2]].x, 
				m_pModel->vpNormals[triangle->NormIdx[2]].y,
				m_pModel->vpNormals[triangle->NormIdx[2]].z);
			glVertex3f(m_pModel->vpVertices[triangle->VertIdx[2]].x,
				m_pModel->vpVertices[triangle->VertIdx[2]].y,
				m_pModel->vpVertices[triangle->VertIdx[2]].z);
			glEnd();
		}

		group = group->next;
	}
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glPopMatrix();
}

void CAccessObj::DrawGrid(){
	glColor3f(1, 1, 0.3);
	glBegin ( GL_LINES );	
	float TCELL = CellSize;
	for (float z = m_vMin.z; z <= m_vMax.z; z+=TCELL ) {
		for (float y = m_vMin.y; y <= m_vMax.y; y+=TCELL ) {
			glVertex3f ( m_vMin.x, y, z );	glVertex3f ( m_vMax.x, y, z );
		}
	}
	for (float z=m_vMin.z; z <= m_vMax.z; z+=TCELL) {
		for (float x=m_vMin.x; x <= m_vMax.x; x+=TCELL) {
			glVertex3f ( x, m_vMin.y, z );	glVertex3f ( x, m_vMax.y, z );
		}
	}
	for (float y=m_vMin.y; y <= m_vMax.y; y+=TCELL ) {
		for (float x=m_vMin.x; x <= m_vMax.x; x+=TCELL ) {
			glVertex3f ( x, y, m_vMin.z );	glVertex3f ( x, y, m_vMax.z );
		}
	}
	glEnd ();
}

//////////////////////////////////////////////////////////////////////
// OpenGLList: Generates and returns a display list for the model using
// the mode specified.
//
// model    - initialized COBJmodel structure
// mode     - a bitwise OR of values describing what is to be rendered.
//            OBJ_NONE     -  render with only vertices
//            OBJ_FLAT     -  render with facet normals
//            OBJ_SMOOTH   -  render with vertex normals
//            OBJ_TEXTURE  -  render with texture coords
//            OBJ_COLOR    -  render with colors (color material)
//            OBJ_MATERIAL -  render with materials
//            OBJ_COLOR and OBJ_MATERIAL should not both be specified.  
// OBJ_FLAT and OBJ_SMOOTH should not both be specified.
//////////////////////////////////////////////////////////////////////
unsigned int CAccessObj::OpenGLList()
{
	unsigned int list;
	
	list = glGenLists(1);
	glNewList(list, GL_COMPILE);
	Draw();
	glEndList();
	
	return list;
}

void CAccessObj::Boundingbox(Vector3DF &vMax, Vector3DF &vMin)
{
	vMax = m_vMax;
	vMin = m_vMin;
}

void CAccessObj::CalcBoundingBox()
{
	unsigned int i;

	m_vMax = m_vMin = m_pModel->vpVertices[1];
	for (i = 1; i <= m_pModel->nVertices; i++)
	//m_vMax = m_vMin = m_pModel->vpVertices[0];
	//for (i = 0; i < m_pModel->nVertices; i++)
	{
		if (m_vMax.x < m_pModel->vpVertices[i].x)
			m_vMax.x = m_pModel->vpVertices[i].x;
		if (m_vMin.x > m_pModel->vpVertices[i].x)
			m_vMin.x = m_pModel->vpVertices[i].x;
		if (m_vMax.y < m_pModel->vpVertices[i].y)
			m_vMax.y = m_pModel->vpVertices[i].y;
		if (m_vMin.y > m_pModel->vpVertices[i].y)
			m_vMin.y = m_pModel->vpVertices[i].y;
		
		if (m_vMax.z < m_pModel->vpVertices[i].z)
			m_vMax.z = m_pModel->vpVertices[i].z;
		if (m_vMin.z > m_pModel->vpVertices[i].z)
			m_vMin.z = m_pModel->vpVertices[i].z;
	}
}

void CAccessObj::FirstGroup()
{
	m_pCurGroup = m_pModel->pGroups;
}

void CAccessObj::NextGroup()
{
	if (m_pCurGroup!=NULL)
		m_pCurGroup = m_pCurGroup->next;
}

unsigned int * CAccessObj::GetGroupTriPIdx()
{
	return m_pCurGroup->pTriangles;
}

int CAccessObj::GetGroupNumTri()
{
	return m_pCurGroup->nTriangles;
}

void CAccessObj::WritePoly(char *filename)
{
    // calculate the face normals
    FacetNormals();

	VertexNormals( 90.0 );
    Unitize();

    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
    // calculate normal of each vertex
    Vector3DF *pVertNormal;

    // create the normals of each vertex
    pVertNormal = new Vector3DF [ m_pModel->nVertices + 1 ];
	//pVertNormal = new Vector3DF [ m_pModel->nVertices ];

    // for each triangle
	//将每个面在此顶点的法向都加了起来.
    for( unsigned int iTri = 0; iTri < m_pModel->nTriangles; iTri++ )
    {
        // get the triangle
        COBJtriangle *pTriangle = &( T( iTri ) );

        // get the normal of the triangle
        Vector3DF pNormal = m_pModel->vpFacetNorms[ pTriangle->findex ];

        // add the normal to each vertex of the triangle
        pVertNormal[ pTriangle->VertIdx[0] ] = pVertNormal[ pTriangle->VertIdx[0] ] + pNormal;
        pVertNormal[ pTriangle->VertIdx[1] ] = pVertNormal[ pTriangle->VertIdx[1] ] + pNormal;
        pVertNormal[ pTriangle->VertIdx[2] ] = pVertNormal[ pTriangle->VertIdx[2] ] + pNormal;
    }

    // average of each vertex normal
    for( unsigned int iVert = 1; iVert <= m_pModel->nVertices; iVert++)
	//for( unsigned int iVert = 0; iVert < m_pModel->nVertices; iVert++)
    {
        pVertNormal[ iVert ].Normalize();
    }


    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*

    FILE *fp = fopen( filename, "wt" );
    if( fp == NULL ) 
        assert( 0 );

    // count the number of triangles ( faces )
    COBJgroup *group = m_pModel->pGroups;
    int numFace = 0;
    while( group )
    {
        numFace += group->nTriangles;
        group = group -> next;
    }
    
    // write the header
    fprintf( fp, "vertices: %d\n", m_pModel->nVertices );
    fprintf( fp, "faces: %d\n", numFace);

    // write vertices and normals(default)
    for(unsigned int i=1; i<= m_pModel->nVertices; i++)
	//for(unsigned int i=0; i< m_pModel->nVertices; i++)
    {
        fprintf( fp, "v %f %f %f %f %f %f\n", 
            m_pModel->vpVertices[i].x,
            m_pModel->vpVertices[i].y,
            m_pModel->vpVertices[i].z,
            pVertNormal[i].x,
            pVertNormal[i].y,
            pVertNormal[i].z );
    }

    delete [] pVertNormal;
	pVertNormal = NULL;

    // write faces 
    group = m_pModel->pGroups;
    while( group )
    {
        for(unsigned int i = 0; i<group->nTriangles; i++)
        {
            fprintf( fp, "f %d %d %d\n", 
                T( group->pTriangles[i] ).VertIdx[0],
                T( group->pTriangles[i] ).VertIdx[1],
                T( group->pTriangles[i] ).VertIdx[2]  );
        }
        group = group -> next;
    }

	//也就是说，文件中只有点、点的法向及面，所有group已合并成一个.
    fclose( fp );
}

void CAccessObj::CorrectMesh(float dis){
	int nTri = m_pModel->nTriangles;
	int nVer = m_pModel->nVertices;

	COBJtriangle* pTri = m_pModel->pTriangles;
	Vector3DF *Vert = m_pModel->vpVertices;
	Vector3DF gc[3];
	Vector3DF pos[3];
	int       TriVer[3];
	Vector3DF a, b, c;


	for(int i = 0; i < m_pModel->nTriangles; i++){
		//一个三角形的三个顶点
		for(int j = 0; j < 3; j++){
			TriVer[j] = pTri->VertIdx[j];
			pos[j] = *(Vert + TriVer[j]);

		}

		if((pos[0]-pos[1]).Length() > 4*dis){
			a = (pos[0]+pos[1])*(float)0.5;
			nVer++;
			nTri++;	
		}
	}

}

void  CAccessObj::initTriGrid()
{
	int x, y, z;

	//CalcBoundingBox();

	/*if(m_vMin.x < 0){
		m_vMin.x = (int)(m_vMin.x-0.5);
	}else m_vMin.x = (int)(m_vMin.x+0.5);

	if(m_vMin.y < 0){
		m_vMin.y = (int)(m_vMin.y-0.5);
	}else m_vMin.y = (int)(m_vMin.y+0.5);

	if(m_vMin.z < 0){
		m_vMin.z = (int)(m_vMin.z-0.5);
	}else m_vMin.z = (int)(m_vMin.z+0.5);*/

	//m_vMin.x = (int)(m_vMin.x-1);
	//m_vMin.y = (int)(m_vMin.y-1);
	//m_vMin.z = (int)(m_vMin.z-1);
	

	x = ceil((m_vMax.x-m_vMin.x)/CellSize);
	y = ceil((m_vMax.y-m_vMin.y)/CellSize);
	z = ceil((m_vMax.z-m_vMin.z)/CellSize);

	m_vMax.x = m_vMin.x + x*CellSize;
	m_vMax.y = m_vMin.y + y*CellSize;
	m_vMax.z = m_vMin.z + z*CellSize;

	m_Res_TreGrid.x = x;
	m_Res_TreGrid.y = y;
	m_Res_TreGrid.z = z;

	nCell = x*y*z;

	GridCell = (int *)malloc(sizeof(int) * (m_pModel->nVertices+1));
	GridTri = (node **)malloc(sizeof(node *) * (nCell+1));

	//第n个粒子所在的网格
	memset(GridCell, 0, sizeof(GridCell));
	//memset(GridTri, NULL, sizeof(GridTri));
	for(int i = 0; i <= nCell; i++)	GridTri[i] = NULL;

	int cell = 0;
	int m_GridSrch = 3;
	m_GridAdjCnt = m_GridSrch*m_GridSrch*m_GridSrch;
	for (int y=0; y < m_GridSrch; y++ ) 
		for (int z=0; z < m_GridSrch; z++ ) 
			for (int x=0; x < m_GridSrch; x++ ) 
				m_GridAdj[cell++] = ( y*m_Res_TreGrid.z + z )*m_Res_TreGrid.x +  x ;			// -1 compensates for ndx 0=empty
	nadj = (m_Res_TreGrid.z + 1)*m_Res_TreGrid.x + 1;

	//高斯插值
	Weight[0] = 0.225;
	Weight[1] = Weight[2] = Weight[3] = 0.13239417;
	Weight[4] = Weight[5] = Weight[6] = 0.12593917;
	BarycCoord[0] = Vector3DF(0.33333333, 0.33333333, 0.33333333);
	float a = 0.05971587, b = 0.47014206, c = 0.79742699, d = 0.10128651;
	BarycCoord[1] = Vector3DF(a, b, b);
	BarycCoord[2] = Vector3DF(b, a, b);
	BarycCoord[3] = Vector3DF(b, b, a);
	BarycCoord[4] = Vector3DF(c, d, d);
	BarycCoord[5] = Vector3DF(d, c, d);
	BarycCoord[6] = Vector3DF(d, d, c);
	//红杉密度为440kg/m3， 设粒子半径同样为0.02，单位面积的粒子的质量为：m=V*D = 0.2*440 = 88kg
	//为了修正水的密度，所以应该用水的密度作为刚体密度
	WUnitMass = 88;
	WoodDensity = 1.0/1000;

	findedTri = (int *)malloc(sizeof(int)*(m_pModel->nTriangles+1));
	memset(findedTri, 0, sizeof(findedTri));

}
void  CAccessObj::InsertGrid()
{
	int nTri = m_pModel->nTriangles;
	COBJtriangle* pTri = m_pModel->pTriangles;
	Vector3DF *Vert = m_pModel->vpVertices;
	Vector3DF gc[3];
	Vector3DF pos[3];

	int TriVer[3];
	int gs[3] = {-1, -1, -1};

	for(int i = 0; i < nTri; i++){
		//一个三角形的三个顶点
		int flag = 1;
		for(int j = 0; j < 3; j++){
			TriVer[j] = pTri->VertIdx[j];
			pos[j] = *(Vert + TriVer[j]);
			//gc[j]代表j顶点所在的网格的x,y,z编号
			gs[j] = ComputPos(pos[j], gc[j]);
			//舍弃网格以外的三角面片
			if(gc[j].x < 0 || gc[j].y < 0 || gc[j].z < 0){ flag = 0; break; }
			if(gc[j].x >= m_Res_TreGrid.x || gc[j].y >= m_Res_TreGrid.y || gc[j].z >= m_Res_TreGrid.z){ flag = 0; break; }

			GridCell[TriVer[j]] = gs[j];	
			if(gs[j] != gs[(j+1)%3] && gs[j] != gs[(j+2)%3]){
				node *pNode = (node *)malloc(sizeof(node));
				pNode->next = NULL;

				pNode->TriNum = i;
				pNode->next = GridTri[gs[j]];
				GridTri[gs[j]] = pNode;
			}
		
		}


		//测量三角形贯穿的其他网格
		if(flag == 1){
			fillGrid(gc[0], gc[1], pos[0], pos[1], i);
			fillGrid(gc[0], gc[2], pos[0], pos[2], i);
			fillGrid(gc[1], gc[2], pos[1], pos[2], i);
		}


		pTri++;

	}

}

int  CAccessObj::ComputPos(Vector3DF pos, Vector3DF &Res)
{
	 int x, y, z;

	x = (int)( (pos.x - m_vMin.x) * gridDelta.x);			// Cell in which particle is located
	y = (int)( (pos.y - m_vMin.y) * gridDelta.y);
	z = (int)( (pos.z - m_vMin.z) * gridDelta.z);

	//碰到碰撞盒最大值的边界的情况
	if(x == m_Res_TreGrid.x) x -= 1;
	if(y == m_Res_TreGrid.y) y -= 1;
	if(z == m_Res_TreGrid.z) z -= 1;

	Res.x = x;
	Res.y = y;
	Res.z = z;

	return (int)( (y*m_Res_TreGrid.z + z)*m_Res_TreGrid.x + x);		
	
}

int  CAccessObj::ComputPos(float *ipos, float *res)
{
	Vector3DF pos;
	Vector3DF Res;
	int x, y, z;
	int gs;

	pos.x = ipos[0];
	pos.y = ipos[1];
	pos.z = ipos[2];

	gs = ComputPos(pos, Res);	
	res[0] = Res.x;
	res[1] = Res.y;
	res[2] = Res.z;

	return gs;
}


float CAccessObj::dist_Point2(Vector3DF a, Vector3DF b){
	float x, y, z;
	x = a.x-b.x;
	y = a.y-b.y;
	z = a.z-b.z;
	return x*x + y*y + z*z;
}
Vector3DF CAccessObj::closest_Tri(int TriNum, Vector3DF p){
	float R[3];
	Vector3DF r;
	int nTri = m_pModel->nTriangles;
	COBJtriangle* pTri = m_pModel->pTriangles+TriNum;
	Vector3DF *Vert = m_pModel->vpVertices;
	Vector3DF a, b, c;

	a = *(Vert + pTri->VertIdx[0]);
	b = *(Vert + pTri->VertIdx[1]);
	c = *(Vert + pTri->VertIdx[2]);

	Vector3DF ab = b-a;
	Vector3DF ac = c-a;
	Vector3DF bc = c-b;

	//s, t, u分别对应ab, ac, bc 的重心坐标, sde是一条边的另一个顶点的方向
	//边ab
	float snom = (p-a).Dot(ab);
	float sdenom = (p-b).Dot(a-b);
	//边ac
	float tnom = (p-a).Dot(ac);
	float tdenom = (p-c).Dot(a-c);

	if(snom <= 0.0 && tnom <= 0.0) { return a; }
	//边bc
	float unom = (p-b).Dot(bc);
	float udenom = (p-c).Dot(b-c);

	if(sdenom <= 0.0 && unom <= 0.0) { return b; }
	if(tdenom <= 0.0 && udenom <= 0.0) { return c; }

	//开始计算edge Voronoi region
	Vector3DF n = (b-a);
	n.Cross(c-a);
	float vc = n.Dot((a-p).Cross(b-p));
	if(vc < 0.0 && snom >= 0.0 && sdenom >= 0.0) { return  a + ab*(snom/(snom+sdenom));	}

	float va = n.Dot((b-p).Cross(c-p));
	if(va < 0.0 && unom >= 0.0 && udenom >= 0.0) { return  b + bc*(unom/(unom+udenom)); }
	
	float vb = n.Dot((c-p).Cross(a-p));
	if(vb < 0.0 && tnom >= 0.0 && tdenom >= 0.0) { return  a + ac*(tnom/(tnom+tdenom)); }


	//否则，最近点在三角形内部
	float u = va/(va+vb+vc);
	float v = vb/(va+vb+vc);
	float w = vc/(va+vb+vc);
	if(abs(u+v+w-1) > 0.1){
		return a;
	}

	Vector3DF inside = a*u + b*v + c*w;

	return inside;
}
float* CAccessObj::closest_Tri(int TriNum, double* Pos){
	Vector3DF p = Vector3DF(Pos[0], Pos[1], Pos[2]);
	Vector3DF R = closest_Tri(TriNum, p);
	float r[3];
	r[0] = R.x; r[1] = R.y; r[2] = R.z;
	return r;
}

//计算包含这个三角形的一部分但是没有顶点的cell
//p1和p2是一个直线的两个端点所在的网格的编号坐标，nTri是三角形的编码，从0到n
//若这个三角形的某一个边经过一个网格，就把给这个网格的链表增加一个表示这个三角形的node
void  CAccessObj::fillGrid(Vector3DF p1, Vector3DF p2, Vector3DF pos1, Vector3DF pos2, int nTri){
	Vector3DF pmin, pmax;
	int gs;

	//if((p1.x == p2.x) && (p1.y == p2.y) && (p1.z == p2.z))return;
	if(abs(p1.x-p2.x)< 2 && abs(p1.y-p2.y)<2 && abs(p1.z-p2.z)<2){
		return;
	}

	//两个端点形成一个长方体的绘制区域
	pmin.x = min(p1.x, p2.x);    pmin.y = min(p1.y, p2.y);    pmin.z = min(p1.z, p2.z);
	pmax.x = max(p1.x, p2.x);    pmax.y = max(p1.y, p2.y);    pmax.z = max(p1.z, p2.z);   
	for(int x = pmin.x; x <= pmax.x; x++){
		for(int y = pmin.y; y <= pmax.y; y++){
			for(int z = pmin.z; z <= pmax.z; z++){
				
				gs = (y*m_Res_TreGrid.z + z)*m_Res_TreGrid.x + x;

				//若这个三角形已经存在于这个链表里，就不用检测了
				if(GridTri[gs]!=NULL){
					node *p = GridTri[gs];
					int flag = 0;
					while(p){
						if(p->TriNum == nTri) {
							flag = 1;
							break;
						}
						p = p->next;

					}
					if(flag == 1) continue;

				}

				Vector3DF CellPos;
				CellPos.x = x; CellPos.y = y; CellPos.z = z;
				//如果这个边穿过网格CellPos,就把这个三角形加到这个网格的cell里

				if(CohenSutherland(CellPos, pos1, pos2)){
					node *pNode = (node *)malloc(sizeof(node));
					pNode->next = NULL;
					
					pNode->TriNum = nTri;
					pNode->next = GridTri[gs];
					GridTri[gs] = pNode;

				}

			}
		}
	}


}

//p1和p2分别是两个端点的坐标
bool CAccessObj::CohenSutherland(Vector3DF CellPos, Vector3DF p1, Vector3DF p2){
	Vector3DF cmin, cmax;
	int code1[6], code2[6];//
	int num1, num2;
	int mid;

	cmin.x = CellPos.x * CellSize + m_vMin.x;
	cmin.y = CellPos.y * CellSize + m_vMin.y;
	cmin.z = CellPos.z * CellSize + m_vMin.z; 

	cmax.x = cmin.x + CellSize;
	cmax.y = cmin.y + CellSize;
	cmax.z = cmin.z + CellSize;

	memset(code1, 0, sizeof(code1));
	memset(code2, 0, sizeof(code2));

	if(p1.x < cmin.x) code1[0] = 1; else if(p1.x > cmax.x) code1[3] = 1;
	if(p1.y < cmin.y) code1[1] = 1; else if(p1.y > cmax.y) code1[4] = 1;
	if(p1.z < cmin.z) code1[2] = 1; else if(p1.z > cmax.z) code1[5] = 1;

	if(p2.x < cmin.x) code2[0] = 1; else if(p2.x > cmax.x) code2[3] = 1;
	if(p2.y < cmin.y) code2[1] = 1; else if(p2.y > cmax.y) code2[4] = 1;
	if(p2.z < cmin.z) code2[2] = 1; else if(p2.z > cmax.z) code2[5] = 1;

	mid = 1;
	num1 = num2 = 0;
	for(int i = 5; i >= 0; i--){
		num1 += code1[i]*mid;
		num2 += code2[i]*mid;

		mid *= 2;
	}
	int re = num1&num2;
	if(re == 0) return true;

	return false;
}

bool CAccessObj::IntersectSegmentTriangle(Vector3DF p, Vector3DF q, Vector3DF a, Vector3DF b, Vector3DF c, Vector3DF &Cp, float &t){
	float u, v, w;
	float *cp = false;
	Vector3DF ab = b-a;
	Vector3DF ac = c-a;
	Vector3DF qp = p-q;

	Vector3DF n = ab;
	n.Cross(ac);
	Vector3DF ap = p-a;

	float d = qp.Dot(n);
	if(d <= 0.0)return false;

	t = ap.Dot(n);
	if(t < 0 || t > d) return false;

	Vector3DF e = qp;
	e.Cross(ap);
	v = ac.Dot(e);
	if(v < 0 || v > d) return false;
	w = ab.Dot(e) * (-1);
	if(w < 0 || w+v > d) return false;

	float ood = 1.0f/d;
	t *= ood;
	v *= ood;
	w *= ood;
	u = 1-v-w;

	Vector3DF CollideP = a + ab * v + ac * w;
	CollideP = CollideP+ n.Normalize()*(float)0.00001;
	Vector3DF pTest = p - qp*t;
	//if((pTest-CollideP).Length()> 0.01){
	//	Cp = p;
	//}

	Cp = CollideP;
	return true;

}


bool CAccessObj::IntersectTest(Vector3DF prePos, Vector3DF pos, Vector3DF &n, Vector3DF &Cp,  float &t){
	Vector3DF res;
	int gs = ComputPos(pos, res);
	if(gs > nCell || gs < 0)return false;
	if(GridTri[gs] == NULL) return false;
	//网格中的三角形
	node * pnode = GridTri[gs];
	Vector3DF closePos ;
	float tt;


	Vector3DF cn;
	Vector3DF p = prePos;
	Vector3DF q = pos;
	Vector3DF pq = q-p;
	Vector3DF va;

	while(pnode != NULL){
		n = m_pModel->vpFacetNorms[pnode->TriNum];
		if(pq.Dot(n)<0){
			//贯穿检测
			int *TriVer = (int*)m_pModel->pTriangles[pnode->TriNum].VertIdx;
			bool bInter = IntersectSegmentTriangle(p, q, m_pModel->vpVertices[TriVer[0]], m_pModel->vpVertices[TriVer[1]], m_pModel->vpVertices[TriVer[2]], closePos, tt);
			if( bInter){
				t = tt;
				Cp = closePos;
				return true;
			}
		}


		pnode = pnode->next;
	}

	return false;

}


bool CAccessObj::Collide_Tree(Vector3DF pos,  Vector3DF &closeP, float &dis){
	Vector3DF res;
	int gs = ComputPos(pos, res);
	if(gs > nCell || gs < 0)return false;
	if(GridTri[gs] == NULL) return false;

	//网格中的三角形
	node * pnode = GridTri[gs];
	Vector3DF R;
	Vector3DF ClosestP(0, 0, 0);
	int testDis = 3;
	float MAXDIS = testDis * testDis * 3 + 1;
	float min_Dis = MAXDIS;
	float Dis;
	while(pnode != NULL){
		//点pos到三角形的最短距离
		R = closest_Tri(pnode->TriNum, pos);

		//最短距离的平方，当cell大小是10时，最大可以是300
		Dis = dist_Point2(R, pos);
		if(min_Dis >= Dis){
			min_Dis = Dis;
			ClosestP = R;
		}
		pnode = pnode->next;
	}

	if(min_Dis >= MAXDIS){
		return false;
	}

	dis = min_Dis;
	closeP = ClosestP;
	
	return true;
}

Vector3DF CAccessObj::computeRepAdh(Vector3DF pos, Vector3DF v, float h, float scale, float space){

	Vector3DF res;
	Vector3DF force(0, 0, 0);
	int gs = ComputPos(pos, res);
	if(gs > nCell || gs < 0)return force;
	if(GridTri[gs] == NULL) return force;
	//当r为0时取得的力的最大值
	k = 9000000;

	//网格中的三角形
	node * pnode = GridTri[gs];
	Vector3DF R;

	Vector3DF triForce(0, 0, 0);
	float dis2, dis;
	float r2, r;
	float s2 = scale*scale;
	r0 = h;

	float coor_h = h/scale;
	Vector3DF cellPos;
	int planenum = m_Res_TreGrid.x * m_Res_TreGrid.z;	

	for(int i = 0; i < m_pModel->nTriangles; i++) findedTri[i] = 0;

	for (int cell=0; cell < m_GridAdjCnt; cell++) {
		int j = gs - nadj + m_GridAdj[cell] ;
		if(j < 0 || j > nCell)continue;

		int y = j/planenum;
		int z = (j%planenum)/m_Res_TreGrid.x;
		int x = (j%planenum)%(int)m_Res_TreGrid.x;


		//确定需不需要查找隔壁网格中的三角形
		cellPos.x = m_vMin.x + x*CellSize;
		cellPos.y = m_vMin.y + y*CellSize;
		cellPos.z = m_vMin.z + z*CellSize;
		if( (pos.x-coor_h) > (cellPos.x+CellSize) || (pos.x + coor_h) < cellPos.x)continue;
		if( (pos.y-coor_h) > (cellPos.y+CellSize) || (pos.y + coor_h) < cellPos.y)continue;
		if( (pos.z-coor_h) > (cellPos.z+CellSize) || (pos.z + coor_h) < cellPos.z)continue;

		pnode = GridTri [j];
		while(pnode != NULL){
			//通过和三角形的距离判断是不是需要积分的三角形
			R = closest_Tri(pnode->TriNum, pos);
			dis2 = dist_Point2(R, pos);
			dis = sqrt(dis2);
			dis *= scale;
			if(dis > h){
				pnode = pnode->next;
				continue;
			}
		
			//判断这个三角形是不是被计算过了
			int ntri = pnode->TriNum;
			if(*(findedTri + ntri) == 1){
				pnode = pnode->next;
				continue;
			}else{
				*(findedTri + ntri) = 1;
			}
			//取出三个点
			int *ver = (int*)m_pModel->pTriangles[ntri].VertIdx;
			Vector3DF a = m_pModel->vpVertices[ver[0]];
			Vector3DF b = m_pModel->vpVertices[ver[1]];
			Vector3DF c = m_pModel->vpVertices[ver[2]];
			//计算面积
			Vector3DF ab = b-a;
			Vector3DF ac = c-a;
			Vector3DF s = ab;
			s.Cross(ac);
			float A = 0.5*s.Length();

			//得到七个采样点的坐标
			Vector3DF sample[7];
			for(int i = 0; i < 7; i++)
				sample[i] = a*BarycCoord[i].x + b*BarycCoord[i].y + c*BarycCoord[i].z;

			triForce = m_pModel->vpFacetNorms[ntri];
			float forceScale = 0;
			//对三角形的每个点，都计算它对粒子p的受力，包括计算两点的距离和traction
			for(int i = 0; i < 7; i++){
				r2 = dist_Point2(sample[i], pos);
				r = sqrt(r2);
				r *= scale;
				if(r >= h)continue;
				
				//计算T
				float cc = (h-r)*(h-r);
				float c0 = (h-r0)*(h-r0);
				float Tra = cc*(cc-c0);
				Tra *= 1.0/(h*h*r0*(2*h-r0));
				Tra *= k;

				//计算积分然后乘以三角形面积
				forceScale += Weight[i]*Tra;
			}
			triForce = triForce*forceScale*A;
			force = force + triForce;

			pnode = pnode->next;
		}
	}

	return force;
}

Vector3DF CAccessObj::computFric(Vector3DF pos, Vector3DF v,  float h, float scale, float idensity,  float m_LapKern, float  visc, float space){

	Vector3DF res;
	Vector3DF force(0, 0, 0);
	int gs = ComputPos(pos, res);
	if(gs > nCell || gs < 0)return force;
	if(GridTri[gs] == NULL) return force;

	//网格中的三角形
	node * pnode = GridTri[gs];
	Vector3DF R;

	Vector3DF triFricTerm(0, 0, 0);
	Vector3DF triAttrTerm(0 ,0, 0);
	float * tforce = NULL;
	float dis2, dis;
	float r2, r;
	float dx, dy, dz;
	float s2 = scale*scale;
	float batta = 1;
	float jdensity = WoodDensity;	

	float coor_h = h/scale;
	Vector3DF cellPos;
	int planenum = m_Res_TreGrid.x * m_Res_TreGrid.z;	

	for(int i = 0; i < m_pModel->nTriangles; i++) findedTri[i] = 0;

	for (int cell=0; cell < m_GridAdjCnt; cell++) {
		int j = gs - nadj + m_GridAdj[cell] ;
		if(j < 0 || j > nCell)continue;

		int y = j/planenum;
		int z = (j%planenum)/m_Res_TreGrid.x;
		int x = (j%planenum)%(int)m_Res_TreGrid.x;


		//确定需不需要查找隔壁网格中的三角形
		cellPos.x = m_vMin.x + x*CellSize;
		cellPos.y = m_vMin.y + y*CellSize;
		cellPos.z = m_vMin.z + z*CellSize;
		if( (pos.x-coor_h) > (cellPos.x+CellSize) || (pos.x + coor_h) < cellPos.x)continue;
		if( (pos.y-coor_h) > (cellPos.y+CellSize) || (pos.y + coor_h) < cellPos.y)continue;
		if( (pos.z-coor_h) > (cellPos.z+CellSize) || (pos.z + coor_h) < cellPos.z)continue;

		pnode = GridTri [j];
		while(pnode != NULL){
			//通过和三角形的距离判断是不是需要积分的三角形
			R = closest_Tri(pnode->TriNum, pos);
			dis2 = dist_Point2(R, pos);
			dis = sqrt(dis2);
			dis *= scale;
			if(dis > h){
				pnode = pnode->next;
				continue;
			}

			//判断这个三角形是不是被计算过了
			int ntri = pnode->TriNum;
			int finded = *(findedTri + ntri);
			if(finded == 1){
				pnode = pnode->next;
				continue;
			}else{
				*(findedTri + ntri) = 1;
			}
			//取出三个点
			int *ver = (int*)m_pModel->pTriangles[ntri].VertIdx;
			Vector3DF a = m_pModel->vpVertices[ver[0]];
			Vector3DF b = m_pModel->vpVertices[ver[1]];
			Vector3DF c = m_pModel->vpVertices[ver[2]];
			//计算面积
			Vector3DF ab = b-a;
			Vector3DF ac = c-a;
			Vector3DF s = ab;
			s.Cross(ac);
			float A = 0.5*s.Length();
			A *= 1.0/(space*space);

			//得到七个采样点的坐标
			Vector3DF sample[7];
			for(int i = 0; i < 7; i++)
				sample[i] = a*BarycCoord[i].x + b*BarycCoord[i].y + c*BarycCoord[i].z;
			//计算

			triFricTerm = Vector3DF(0, 0, 0);
			triAttrTerm = Vector3DF(0, 0, 0);
			//对三角形的每个点，都计算它对粒子p的粘性力,计算积分然后乘以三角形面积
			for(int i = 0; i < 7; i++){
				r2 = dist_Point2(sample[i], pos);
				r = sqrt(r2);
				r *= scale;
				if(r > h)continue;

				//光滑核
				if(r < h){				
					float cc = (h-r);
					float dterm = cc * idensity * jdensity;
					float vterm = m_LapKern * visc;//粘性力部分
					Vector3DF samfricforce = (v)*vterm*dterm*(-1) * Weight[i];
					triFricTerm += samfricforce;
				}

				//float Fun = /*0.007*/1.0/(pow(h0, (float)3.25));
				//float tmp ;
				//if(r <= h && 2*r > h){
				//	tmp = -4*r*r/h0+6*r-2*h0;
				//	tmp = pow(tmp, (float)0.25);
				//}else{
				//	tmp = 0;
				//}
				//Fun *= tmp;
				//Vector3DF samattrforce = (sample[i]-pos)*batta*Fun*((float)1.0/r) * Weight[i];
				//	
				//triAttrTerm += samattrforce;
	
			}
			force = force + triFricTerm*A + triAttrTerm*A;
			pnode = pnode->next;
		}
	}


	return force;
}