#ifndef __OUR_GL_H__
#define __OUR_GL_H__

#include "tgaimage.h"
#include "geometry.h"
#include <Eigen/Dense>

extern Matrix Model_Matrix;
extern Matrix View_Matrix;
extern Matrix Projection_Matrix;
extern Eigen::Matrix4f Model_Matrix1;
extern Eigen::Matrix4f View_Matrix1;
extern Eigen::Matrix4f Projection_Matrix1;

void set_model(float angle_z, Vec3f transformation, float scale);
void set_view(Vec3f eye, Vec3f center, Vec3f up);
void set_projection(float eye_fov, float aspect, float zNear, float zFar); 

struct IShader {
    virtual ~IShader();
    virtual Vec3f vertex(int iface, int nthvert) = 0;
    virtual float fragment(Vec3f bar, TGAColor& color) = 0;
};

void line(Vec3i p0, Vec3i p1, TGAImage& image, TGAColor color);

void triangle(Vec3f* pts, IShader& shader, TGAImage& image, float* zbuffer);

#endif //__OUR_GL_H__
