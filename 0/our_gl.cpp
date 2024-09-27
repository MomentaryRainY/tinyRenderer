#include <cmath>
#include <limits>
#include <cstdlib>
#include "our_gl.h"


Matrix Model_Matrix;
Matrix View_Matrix;
Matrix Projection_Matrix;
Eigen::Matrix4f Model_Matrix1;
Eigen::Matrix4f View_Matrix1;
Eigen::Matrix4f Projection_Matrix1;

IShader::~IShader() {}

const float PI = 3.1415927;

Matrix rotate_z(float angle) {
	Matrix m = Matrix::identity();
	float theta = angle * PI / 180;
	m[0][0] = m[1][1] = cos(theta);
	m[0][1] = -sin(theta);
	m[1][0] = sin(theta);
	return m;
}

void set_view(Vec3f eye, Vec3f center, Vec3f up) {
	View_Matrix = Matrix::identity();
	Vec3f z = (center - eye).normalize();
	Vec3f x = cross(up, z).normalize();
	Vec3f y = cross(z, x).normalize();
	for (int i = 0; i < 3; i++) {
		View_Matrix[0][i] = x[i];
		View_Matrix[1][i] = y[i];
		View_Matrix[2][i] = z[i];
	}
	View_Matrix[0][3] -= x * eye;
	View_Matrix[1][3] -= y * eye;
	View_Matrix[2][3] -= z * eye;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			View_Matrix1(i, j) = View_Matrix[i][j];
		}
	}
}

void set_projection(float eye_fov, float aspect, float zNear, float zFar) {
	float t = tan((eye_fov * PI / 180) / 2) * fabs(zNear);
	float r = t * aspect;
	float b = -t;
	float l = -r;

	Matrix ndc = Matrix::identity();
	/*ndc[0][0] = 2 / (r - l);
	ndc[1][1] = 2 / (t - b);
	ndc[2][2] = 2 / (zNear - zFar);
	ndc[0][3] = -(r + l) / (r - l);
	ndc[1][3] = -(t + b) / (t - b);
	ndc[2][3] = -(zNear + zFar) / (zNear - zFar);

	Matrix perspective = Matrix::identity();
	perspective[0][0] = zNear;
	perspective[1][1] = zNear;
	perspective[2][2] = zNear + zFar;
	perspective[2][3] = -zNear * zFar;
	perspective[3][2] = -1;
	perspective[3][3] = 0;*/

	Matrix perspective = Matrix::identity();
	float tanfov = 1 / tan((eye_fov * PI / 180) / 2);
	perspective[0][0] = tanfov / aspect;
	perspective[1][1] = -tanfov;
	perspective[2][2] = -(zFar + zNear) / (zFar - zNear);
	perspective[2][3] = -2 * zNear * zFar / (zFar - zNear);
	perspective[3][2] = -1;
	perspective[3][3] = 0; 

	Projection_Matrix = ndc * perspective;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Projection_Matrix1(i, j) = Projection_Matrix[i][j];
		}
	}
}

void set_model(float angle_z, Vec3f transformation, float scale) {
	Matrix trans = Matrix::identity();
	for (int i = 0; i < 3; i++) trans[i][3] = transformation[i];
	Matrix rotation = Matrix::identity();
	for (int i = 0; i < 3; i++) rotation = rotate_z(angle_z);
	Matrix scaleM = Matrix::identity();
	for (int i = 0; i < 3; i++) scaleM[i][i] = scale;
	Model_Matrix = trans * rotation * scaleM;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Model_Matrix1(i, j) = Model_Matrix[i][j];
		}
	}
}

void line(Vec3i p0, Vec3i p1, TGAImage& image, TGAColor color) {
	bool steep = false;
	if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y)) {
		std::swap(p0.x, p0.y);
		std::swap(p1.x, p1.y);
		steep = true;
	}
	if (p0.x > p1.x) {
		std::swap(p0, p1);
	}

	for (int x = p0.x; x <= p1.x; x++) {
		float t = (x - p0.x) / (float)(p1.x - p0.x);
		int y = p0.y * (1. - t) + p1.y * t + .5;
		if (steep) {
			image.set(y, x, color);
		}
		else {
			image.set(x, y, color);
		}
	}
}

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
	Vec3f s[2];
	for (int i = 2; i--; ) {
		s[i][0] = C[i] - A[i];
		s[i][1] = B[i] - A[i];
		s[i][2] = A[i] - P[i];
	}
	Vec3f u = cross(s[0], s[1]);
	if (std::abs(u[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
		return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
	return Vec3f(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void triangle(Vec3f* pts, IShader &shader, TGAImage& image, float* zbuffer) {
	std::ofstream log("log.txt", std::ios::app);
	log << pts[0][2] << " " << pts[1][2] << " " << pts[2][2] << std::endl;
	Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			bboxmin[j] = std::min(bboxmin[j], pts[i][j] );
			bboxmax[j] = std::max(bboxmax[j], pts[i][j] );
		}
	}
	Vec2i P;
	TGAColor color;
	//log << int(color.bgra[0]) << " " << int(color.bgra[1]) << " " << int(color.bgra[2]) << std::endl;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f c = barycentric(pts[0], pts[1], pts[2], Vec3f(P.x, P.y, 1.f));
			float z = pts[0][2] * c.x + pts[1][2] * c.y + pts[2][2] * c.z;
			float frag_depth = std::max(0.f, std::min(255.f, z));
			if (c.x < 0 || c.y < 0 || c.z < 0 || zbuffer[P.x + P.y * 800] > frag_depth ) continue;
			float discard = shader.fragment(c, color);
			//log << "locate " << P.x << " " << P.y << std::endl;
			//log << "intensity " << discard << std::endl;
			if (!false) {
				//log << "buffer " << zbuffer[P.x + P.y * 800] << std::endl;
				zbuffer[P.x + P.y * 800] = frag_depth;
				image.set(P.x, P.y, color);
				//log << "frag " << frag_depth << std::endl;
				//log << int(color.bgra[0]) << " " << int(color.bgra[1]) << " " << int(color.bgra[2]) << std::endl;
				//log << std::endl;
			}
		}
	}
}