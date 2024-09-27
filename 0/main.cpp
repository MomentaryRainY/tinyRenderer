#include <vector>
#include <iostream>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h";


Model* model = NULL;
const int width = 800;
const int height = 800;

Vec3f light_dir = Vec3f(1,1,1).normalize();

/*model matrix parameters*/
float rotateZ = 0;
Vec3f translate(0, 0, 0);
float scale = 1;
/*view matrix parameters*/
Vec3f eye(1, 1, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);
/*projection matrix parameters*/
float eye_fov = 45;
float aspect = width / height;
float zNear = .5f;
float zFar = 3.f;

float f1 = (zFar - zNear) / 2.0;
float f2 = (zFar + zNear) / 2.0;

float* zbuffer = new float[height * width];
float* depthbuffer = new float[height * width];

float M_PI = 3.1415926;

Matrix viewport;

struct ShadowShader : public IShader {
	mat<3, 3, float> varying_tri;

	ShadowShader() : varying_tri() {}

	virtual Vec3f vertex(int iface, int nthvert) {
		Vec4f gl_vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
		gl_vertex = Projection_Matrix * View_Matrix * Model_Matrix * gl_vertex;          // transform it to screen coordinates
		gl_vertex[0] = (gl_vertex[0] / gl_vertex[3] + 1.f) * width / 2.f;
		gl_vertex[1] = (gl_vertex[1] / gl_vertex[3] + 1.f) * height / 2.f;
		gl_vertex[2] = (gl_vertex[2] / gl_vertex[3] - zNear) / (zFar - zNear) * 255;
		//gl_vertex[2] = (gl_vertex[3] / gl_vertex[2] - 1 / zNear) * 255 / (1 / zFar - 1 / zNear);
		varying_tri.set_col(nthvert, Vec3f(gl_vertex[0], gl_vertex[1], gl_vertex[2]));
		return Vec3f(gl_vertex[0], gl_vertex[1], gl_vertex[2]);
	}

	virtual float fragment(Vec3f bar, TGAColor& color) {
		Vec3f p = varying_tri * bar;
		//color = TGAColor(255, 255, 255) * (p.z / 255);
		color = TGAColor(0, 0, 0);
		return p.z;
	}
};

struct GourauShader : public IShader {
	//Vec3f varying_intensity;
	mat<2, 3, float> varying_uv;
	mat<3, 3, float> varying_tri;
	mat<4, 4, float> uniform_M;
	mat<4, 4, float> uniform_MIT;
	mat<4, 4, float> uniform_Mshadow;

	virtual Vec3f vertex(int iface, int nthvert) {
		Vec4f gl_vertex = embed<4>(model->vert(iface, nthvert));
		gl_vertex = Projection_Matrix * View_Matrix * Model_Matrix * gl_vertex;
		gl_vertex[0] = (gl_vertex[0] / gl_vertex[3] + 1.f) * width / 2.f;
		gl_vertex[1] = (gl_vertex[1] / gl_vertex[3] + 1.f) * height / 2.f;
		gl_vertex[2] = (gl_vertex[2] / gl_vertex[3] - zNear) / (zFar - zNear) * 255;
		//varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert) * light_dir);
		varying_uv.set_col(nthvert, model->uv(iface, nthvert));
		varying_tri.set_col(nthvert, Vec3f(gl_vertex[0], gl_vertex[1], gl_vertex[2]));
		return Vec3f(gl_vertex[0], gl_vertex[1], gl_vertex[2]);
	}

	virtual float fragment(Vec3f bar, TGAColor &color) {
		Vec4f sb_p = uniform_Mshadow * embed<4>(varying_tri * bar); // corresponding point in the shadow buffer
		sb_p = sb_p / sb_p[3];
		int idx = int(sb_p[0]) + int(sb_p[1]) * width; // index in the shadowbuffer array
		float shadow = .3 + .7 * (depthbuffer[idx] < sb_p[2] + .5f);
		Vec2f uv = varying_uv * bar;
		Vec3f n = proj<3>(uniform_MIT * embed<4>(model->normal(uv))).normalize();
		Vec3f l = proj<3>(uniform_M * embed<4>(light_dir)).normalize();
		Vec3f r = (n * (n * l * 2.f) - l).normalize();
		float spec = r.z > 0.8 ? pow(r.z, model->specular(uv)) : 0;
		float diff = std::max(0.f, n * l);
		TGAColor c = model->diffuse(uv) * diff;
		color = c;
		for (int i = 0; i < 3; i++) color[i] = std::min<float>(c[i] * shadow * (diff + spec), 255);
		return 1;
	}
};

float max_elevation_angle(float* zbuffer, Vec2f p, Vec2f dir) {
	float maxangle = 0;
	for (float t = 0.; t < 1000.; t += 1.) {
		Vec2f cur = p + dir * t;
		if (cur.x >= width || cur.y >= height || cur.x < 0 || cur.y < 0) return maxangle;

		float distance = (p - cur).norm();
		if (distance < 1.f) continue;
		float elevation = zbuffer[int(cur.x) + int(cur.y) * width] - zbuffer[int(p.x) + int(p.y) * width];
		maxangle = std::max(maxangle, atanf(elevation / distance));
	}
	return maxangle;
}

int main(int argc, char** argv) {

	model = new Model("obj/african_head.txt");
	//model = new Model("obj/cube.txt");

	//std::ofstream log("log.txt", std::ios::trunc);
	Matrix viewport = Matrix::identity();
	viewport[0][0] = viewport[0][3] = width / 2;
	viewport[1][1] = viewport[1][3] = height / 2;
	viewport[2][2] = 255 / (zFar - zNear);
	viewport[2][3] = -zNear * 255 / (zFar - zNear);
	Matrix M;
	{
		TGAImage depth(width, height, TGAImage::RGB);
		for (int i = width * height; i--; depthbuffer[i] = -std::numeric_limits<float>::max());

		set_model(rotateZ, translate, scale);
		set_view(Vec3f(2, 2, 2), center, up);
		set_projection(eye_fov, aspect, zNear, zFar);
		M = viewport * Projection_Matrix * View_Matrix * Model_Matrix;

		ShadowShader sShader;

		for (int i = 0; i < model->nfaces(); i++) {
			Vec3f screen_coords[3];
			for (int j = 2; j >= 0; j--) {
				screen_coords[j] = sShader.vertex(i, j);
			}
			triangle(screen_coords, sShader, depth, depthbuffer);
		}

		depth.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		depth.write_tga_file("depth.tga");
	}

	TGAImage image(width, height, TGAImage::RGB);

	 {
		for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

		set_model(rotateZ, translate, scale);
		set_view(eye, center, up);
		set_projection(eye_fov, aspect, zNear, zFar);

		GourauShader shader;
		Matrix mvp = Projection_Matrix * View_Matrix * Model_Matrix;
		shader.uniform_M = mvp;
		shader.uniform_MIT = mvp.invert_transpose();
		shader.uniform_Mshadow = M * (viewport * mvp).invert();
		Matrix x = (viewport * mvp).invert();


		for (int i = 0; i < model->nfaces(); i++) {
			Vec3f screen_coords[3];
			for (int j = 2; j >= 0; j--) {
				screen_coords[j] = shader.vertex(i, j);
			}
			triangle(screen_coords, shader, image, zbuffer);
		}

	}
	 /*
	{
		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				if (zbuffer[x + y * width] < -1e5) continue;
				float total = 0;
				for (float a = 0; a < M_PI * 2 - 1e-4; a += M_PI / 4)
					total += M_PI - max_elevation_angle(zbuffer, Vec2f(x, y), Vec2f(cos(a), sin(a)));
				total /= (M_PI / 2) * 8;
				total = pow(total, 27.f);
				image.set(x, y, image.get(x, y) + TGAColor(255 * total, 255 * total, 255 * total));
			}
		}
	}*/

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");

	delete[] zbuffer;
	delete[] depthbuffer;
	delete model;
	//log.close();
	return 0;
}