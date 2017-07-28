#if defined(_MSC_VER)
#pragma once
#endif

#ifndef OPENGL_UTILS_H
#define OPENGL_UTILS_H

#include <vector>
#include "types.h"

/************************************************************************/
/*							OpenGL routines                             */
/************************************************************************/

void draw_circle2d(const Vec2r& centre, Real rad, int segs);
void draw_grid2d(const Vec2r& origin, Real dx, int nx, int ny);
void draw_box2d(const Vec2r& origin, Real width, Real height);
void draw_segmentset2d(const std::vector<Vec2r>& vertices, const std::vector<Vec2ui>& edges);
void draw_segmentset2d(const std::vector<Vec2r>& vertices, const std::vector<Vec2i>& edges);
void draw_points2d(const std::vector<Vec2r>& points);
void draw_points2d_color(const std::vector<Vec2r>& points,const std::vector<Vec3r>& color,Real radius);
void draw_polygon2d(const std::vector<Vec2r>& vertices);
void draw_polygon2d(const std::vector<Vec2r>& vertices, const std::vector<int>& order);
void draw_segment2d(const Vec2r& start, const Vec2r& end);
void draw_arrow2d(const Vec2r& start, const Vec2r& end, Real arrow_head_len);
void draw_arrow3d(const Vec3r& start, const Vec3r& end, Real arrow_head_len);
void draw_grid_data2d(Array2f& data, Vec2r origin, Real dx, bool color = false);
void draw_trimesh2d(const std::vector<Vec2r>& vertices, const std::vector<Vec3ui>& tris);    
   
void draw_trimesh3d(const std::vector<Vec3r>& vertices, const std::vector<Vec3ui>& tris);
void draw_trimesh3d(const std::vector<Vec3r>& vertices, const std::vector<Vec3ui>& tris, const std::vector<Vec3r>& normals);
void draw_box3d(const Vec3r& dimensions);
void draw_line_box3d(const Vec3r& dimensions);
void draw_line_box3d(const Vec3r& left, const Vec3r& right);
void draw_densities(const Vec2r& position, Real dx,Real p1,Real p2,Real p3,Real p4);
void draw_points_cloud(const std::vector<Vec2r>& points,Real radius);

void drawBitmapString( const char *string);

#endif