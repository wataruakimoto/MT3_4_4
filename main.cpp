#include <Novice.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>
#include <imgui.h>

struct Vector3 {
	float x;
	float y;
	float z;
};

struct Matrix4x4 {
	float m[4][4];
};

struct Plane {
	Vector3 normal; // 法線
	float distance; // 距離
};

struct Sphere {
	Vector3 center; // 中心点
	float radius; // 半径
};

struct Ball {
	Vector3 position;     // ボールの位置
	Vector3 velocity;     // ボールの速度
	Vector3 acceleration; // ボールの加速度
	float mass;           // ボールの質量
	float radius;         // ボールの半径
	unsigned int color;   // ボールの色
};

// 加算
Vector3 Add(const Vector3& v1, const Vector3& v2) {

	Vector3 resultAdd;

	resultAdd.x = v1.x + v2.x;
	resultAdd.y = v1.y + v2.y;
	resultAdd.z = v1.z + v2.z;

	return resultAdd;
}

Vector3 operator+(const Vector3& v1, const Vector3& v2) { return Add(v1, v2); }

// 減算
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {

	Vector3 resultSubtract;

	resultSubtract.x = v1.x - v2.x;
	resultSubtract.y = v1.y - v2.y;
	resultSubtract.z = v1.z - v2.z;

	return resultSubtract;
}

Vector3 operator-(const Vector3& v1, const Vector3& v2) { return Subtract(v1, v2); }

// 積
Vector3 Multiply(float scalar, const Vector3& v) {

	Vector3 resultMultiply;

	resultMultiply.x = v.x * scalar;
	resultMultiply.y = v.y * scalar;
	resultMultiply.z = v.z * scalar;

	return resultMultiply;
}

Vector3 operator*(float s, const Vector3& v) { return Multiply(s, v); }
Vector3 operator*(const Vector3& v, float s) { return s * v; }

//内積
float Dot(const Vector3& v1, const Vector3& v2) {

	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

//長さ（ノルム）
float Length(const Vector3& v) {

	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

//正規化
Vector3 Normalize(const Vector3& v) {

	Vector3 result;

	result.x = v.x / sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	result.y = v.y / sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	result.z = v.z / sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);

	return result;
}

// クロス積
Vector3 Cross(const Vector3& v1, const Vector3& v2) {

	Vector3 resultCross = {};

	resultCross.x = v1.y * v2.z - v1.z * v2.y;
	resultCross.y = v1.z * v2.x - v1.x * v2.z;
	resultCross.z = v1.x * v2.y - v1.y * v2.x;

	return resultCross;
}

// 座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {

	// w = 1 がデカルト座標系であるので(x,y,z,1)のベクトルとしてmatrixとの積をとる
	Vector3 resultTransform = {};

	resultTransform.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
	resultTransform.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
	resultTransform.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];

	// ベクトルに対して基本的な操作を行う行列でwが0になることはありえない
	assert(w != 0.0f);

	// w = 1 がデカルト座標系であるので、w除算することで同次座標をデカルト座標に戻す
	resultTransform.x /= w;
	resultTransform.y /= w;
	resultTransform.z /= w;

	return resultTransform;
}

// 行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {

	Matrix4x4 resultMultiply = {};

	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			resultMultiply.m[y][x] = m1.m[y][0] * m2.m[0][x] + m1.m[y][1] * m2.m[1][x] + m1.m[y][2] * m2.m[2][x] + m1.m[y][3] * m2.m[3][x];
		}
	}

	return resultMultiply;
}

// 逆行列
Matrix4x4 Inverse(const Matrix4x4& m) {

	float determinant =
		m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]
		- m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]
		- m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]
		+ m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]
		+ m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]
		- m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]
		- m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]
		+ m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];

	Matrix4x4 resultInverse = {};

	resultInverse.m[0][0] = (m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[1][3] * m.m[2][1] * m.m[3][2]
		- m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2]) / determinant;

	resultInverse.m[0][1] = (-m.m[0][1] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[2][1] * m.m[3][2]
		+ m.m[0][3] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[2][3] * m.m[3][2]) / determinant;

	resultInverse.m[0][2] = (m.m[0][1] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[3][2]
		- m.m[0][3] * m.m[1][2] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[3][2]) / determinant;

	resultInverse.m[0][3] = (-m.m[0][1] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][1] - m.m[0][3] * m.m[1][1] * m.m[2][2]
		+ m.m[0][3] * m.m[1][2] * m.m[2][1] + m.m[0][2] * m.m[1][1] * m.m[2][3] + m.m[0][1] * m.m[1][3] * m.m[2][2]) / determinant;


	resultInverse.m[1][0] = (-m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[1][3] * m.m[2][0] * m.m[3][2]
		+ m.m[1][3] * m.m[2][2] * m.m[3][0] + m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[1][0] * m.m[2][3] * m.m[3][2]) / determinant;

	resultInverse.m[1][1] = (m.m[0][0] * m.m[2][2] * m.m[3][3] + m.m[0][2] * m.m[2][3] * m.m[3][0] + m.m[0][3] * m.m[2][0] * m.m[3][2]
		- m.m[0][3] * m.m[2][2] * m.m[3][0] - m.m[0][2] * m.m[2][0] * m.m[3][3] - m.m[0][0] * m.m[2][3] * m.m[3][2]) / determinant;

	resultInverse.m[1][2] = (-m.m[0][0] * m.m[1][2] * m.m[3][3] - m.m[0][2] * m.m[1][3] * m.m[3][0] - m.m[0][3] * m.m[1][0] * m.m[3][2]
		+ m.m[0][3] * m.m[1][2] * m.m[3][0] + m.m[0][2] * m.m[1][0] * m.m[3][3] + m.m[0][0] * m.m[1][3] * m.m[3][2]) / determinant;

	resultInverse.m[1][3] = (m.m[0][0] * m.m[1][2] * m.m[2][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] + m.m[0][3] * m.m[1][0] * m.m[2][2]
		- m.m[0][3] * m.m[1][2] * m.m[2][0] - m.m[0][2] * m.m[1][0] * m.m[2][3] - m.m[0][0] * m.m[1][3] * m.m[2][2]) / determinant;


	resultInverse.m[2][0] = (m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[1][3] * m.m[2][0] * m.m[3][1]
		- m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1]) / determinant;

	resultInverse.m[2][1] = (-m.m[0][0] * m.m[2][1] * m.m[3][3] - m.m[0][1] * m.m[2][3] * m.m[3][0] - m.m[0][3] * m.m[2][0] * m.m[3][1]
		+ m.m[0][3] * m.m[2][1] * m.m[3][0] + m.m[0][1] * m.m[2][0] * m.m[3][3] + m.m[0][0] * m.m[2][3] * m.m[3][1]) / determinant;

	resultInverse.m[2][2] = (m.m[0][0] * m.m[1][1] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] + m.m[0][3] * m.m[1][0] * m.m[3][1]
		- m.m[0][3] * m.m[1][1] * m.m[3][0] - m.m[0][1] * m.m[1][0] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][1]) / determinant;

	resultInverse.m[2][3] = (-m.m[0][0] * m.m[1][1] * m.m[2][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] - m.m[0][3] * m.m[1][0] * m.m[2][1]
		+ m.m[0][3] * m.m[1][1] * m.m[2][0] + m.m[0][1] * m.m[1][0] * m.m[2][3] + m.m[0][0] * m.m[1][3] * m.m[2][1]) / determinant;


	resultInverse.m[3][0] = (-m.m[1][0] * m.m[2][1] * m.m[3][2] - m.m[1][1] * m.m[2][2] * m.m[3][0] - m.m[1][2] * m.m[2][0] * m.m[3][1]
		+ m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[1][1] * m.m[2][0] * m.m[3][2] + m.m[1][0] * m.m[2][2] * m.m[3][1]) / determinant;

	resultInverse.m[3][1] = (m.m[0][0] * m.m[2][1] * m.m[3][2] + m.m[0][1] * m.m[2][2] * m.m[3][0] + m.m[0][2] * m.m[2][0] * m.m[3][1]
		- m.m[0][2] * m.m[2][1] * m.m[3][0] - m.m[0][1] * m.m[2][0] * m.m[3][2] - m.m[0][0] * m.m[2][2] * m.m[3][1]) / determinant;

	resultInverse.m[3][2] = (-m.m[0][0] * m.m[1][1] * m.m[3][2] - m.m[0][1] * m.m[1][2] * m.m[3][0] - m.m[0][2] * m.m[1][0] * m.m[3][1]
		+ m.m[0][2] * m.m[1][1] * m.m[3][0] + m.m[0][1] * m.m[1][0] * m.m[3][2] + m.m[0][0] * m.m[1][2] * m.m[3][1]) / determinant;

	resultInverse.m[3][3] = (m.m[0][0] * m.m[1][1] * m.m[2][2] + m.m[0][1] * m.m[1][2] * m.m[2][0] + m.m[0][2] * m.m[1][0] * m.m[2][1]
		- m.m[0][2] * m.m[1][1] * m.m[2][0] - m.m[0][1] * m.m[1][0] * m.m[2][2] - m.m[0][0] * m.m[1][2] * m.m[2][1]) / determinant;

	return resultInverse;
}

// 3次元回転行列
Matrix4x4 MakeRotateMatrix(Vector3 radian) {

	Matrix4x4 rotateX = { 0.0f };

	rotateX.m[0][0] = 1.0f;
	rotateX.m[1][1] = cosf(radian.x);
	rotateX.m[1][2] = sinf(radian.x);
	rotateX.m[2][1] = -sinf(radian.x);
	rotateX.m[2][2] = cosf(radian.x);
	rotateX.m[3][3] = 1.0f;

	Matrix4x4 rotateY = { 0.0f };

	rotateY.m[0][0] = cosf(radian.y);
	rotateY.m[0][2] = -sinf(radian.y);
	rotateY.m[1][1] = 1.0f;
	rotateY.m[2][0] = sinf(radian.y);
	rotateY.m[2][2] = cosf(radian.y);
	rotateY.m[3][3] = 1.0f;

	Matrix4x4 rotateZ = { 0.0f };

	rotateZ.m[0][0] = cosf(radian.z);
	rotateZ.m[0][1] = sinf(radian.z);
	rotateZ.m[1][0] = -sinf(radian.z);
	rotateZ.m[1][1] = cosf(radian.z);
	rotateZ.m[2][2] = 1.0f;
	rotateZ.m[3][3] = 1.0f;

	Matrix4x4 resultRotate = { 0.0f };

	resultRotate = Multiply(rotateX, Multiply(rotateY, rotateZ));

	return resultRotate;
}

// 3次元アフィン変換行列
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {

	Matrix4x4 rotateMatrix = { 0.0f };

	rotateMatrix = MakeRotateMatrix(rotate);

	Matrix4x4 resultAffine = { 0.0f };

	resultAffine.m[0][0] = scale.x * rotateMatrix.m[0][0];
	resultAffine.m[0][1] = scale.x * rotateMatrix.m[0][1];
	resultAffine.m[0][2] = scale.x * rotateMatrix.m[0][2];
	resultAffine.m[1][0] = scale.y * rotateMatrix.m[1][0];
	resultAffine.m[1][1] = scale.y * rotateMatrix.m[1][1];
	resultAffine.m[1][2] = scale.y * rotateMatrix.m[1][2];
	resultAffine.m[2][0] = scale.z * rotateMatrix.m[2][0];
	resultAffine.m[2][1] = scale.z * rotateMatrix.m[2][1];
	resultAffine.m[2][2] = scale.z * rotateMatrix.m[2][2];

	resultAffine.m[3][0] = translate.x;
	resultAffine.m[3][1] = translate.y;
	resultAffine.m[3][2] = translate.z;
	resultAffine.m[3][3] = 1.0f;

	return resultAffine;
}

// 透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {

	Matrix4x4 resultPerspectiveFov = {};

	resultPerspectiveFov.m[0][0] = (1 / aspectRatio) * (1 / tanf(fovY / 2));
	resultPerspectiveFov.m[1][1] = 1 / tanf(fovY / 2);
	resultPerspectiveFov.m[2][2] = farClip / (farClip - nearClip);
	resultPerspectiveFov.m[2][3] = 1.0f;
	resultPerspectiveFov.m[3][2] = -nearClip * farClip / (farClip - nearClip);

	return resultPerspectiveFov;
}

// ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {

	Matrix4x4 resultViewport = {};

	resultViewport.m[0][0] = width / 2;
	resultViewport.m[1][1] = -height / 2;
	resultViewport.m[2][2] = maxDepth - minDepth;
	resultViewport.m[3][0] = left + (width / 2);
	resultViewport.m[3][1] = top + (height / 2);
	resultViewport.m[3][2] = minDepth;
	resultViewport.m[3][3] = 1.0f;

	return resultViewport;
}

// 球の描画
void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {

	const uint32_t kSubdivision = 16; // 分割数
	const float kLatEvery = (float)M_PI / kSubdivision; // 経度分割1つ分の角度
	const float kLonEvery = (2 * (float)M_PI) / kSubdivision; // 緯度分割1つ分の角度

	// 緯度の方向に分割 -Π/2 ~ Π/2
	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {

		float lat = -(float)M_PI / 2.0f + latIndex * kLatEvery; // 現在の緯度

		// 経度の方向に分割 0 ~ 2Π
		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {

			float lon = lonIndex * kLonEvery; // 現在の経度

			// World座標系でのa、b、cを求める
			Vector3 a, b, c;

			a = {
				sphere.center.x + sphere.radius * cosf(lat) * cosf(lon),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * cosf(lat) * sinf(lon)
			};

			b = {
				sphere.center.x + sphere.radius * cosf(lat + kLatEvery) * cosf(lon),
				sphere.center.y + sphere.radius * sinf(lat + kLatEvery),
				sphere.center.z + sphere.radius * cosf(lat + kLatEvery) * sinf(lon)
			};

			c = {
				sphere.center.x + sphere.radius * cosf(lat) * cosf(lon + kLonEvery),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * cosf(lat) * sinf(lon + kLonEvery)
			};

			// a、b、cをScreen座標系まで変換
			a = Transform(a, Multiply(viewProjectionMatrix, viewportMatrix));
			b = Transform(b, Multiply(viewProjectionMatrix, viewportMatrix));
			c = Transform(c, Multiply(viewProjectionMatrix, viewportMatrix));

			// ab、bcで線を引く
			Novice::DrawLine(int(a.x), int(a.y), int(b.x), int(b.y), color);
			Novice::DrawLine(int(a.x), int(a.y), int(c.x), int(c.y), color);
		}
	}
}

// 矩形の4頂点を求める
Vector3 Perpendicular(const Vector3& vector) {

	if (vector.x != 0.0f || vector.y != 0.0f) {

		return { -vector.y,vector.x,0.0f };
	}
	return{ 0.0f,-vector.z,vector.y };
}

// 平面の描画
void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {

	// 1:中心点を求める
	Vector3 center = Multiply(plane.distance, plane.normal);

	Vector3 perpendiculars[4];

	// 2:法線と垂直なベクトルを1つ求める
	perpendiculars[0] = Normalize(Perpendicular(plane.normal));

	// 3:2の逆ベクトルを求める
	perpendiculars[1] = { -perpendiculars[0].x,-perpendiculars[0].y,-perpendiculars[0].z };

	// 4:2と法線とのクロス積を求める
	perpendiculars[2] = Cross(plane.normal, perpendiculars[0]);

	// 5:4の逆ベクトルを求める
	perpendiculars[3] = { -perpendiculars[2].x,-perpendiculars[2].y,-perpendiculars[2].z };

	// 6:2~5のベクトルを中心点にそれぞれ定数倍して足すと4頂点が出来上がる
	Vector3 points[4];
	for (int32_t index = 0; index < 4; ++index) {
		Vector3 extend = Multiply(2.0f, perpendiculars[index]);
		Vector3 point = Add(center, extend);
		points[index] = Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
	}

	// pointsをそれぞれ結んでDrawLineで矩形を描画する。
	Novice::DrawLine((int)points[0].x, (int)points[0].y, (int)points[2].x, (int)points[2].y, color);
	Novice::DrawLine((int)points[1].x, (int)points[1].y, (int)points[2].x, (int)points[2].y, color);
	Novice::DrawLine((int)points[1].x, (int)points[1].y, (int)points[3].x, (int)points[3].y, color);
	Novice::DrawLine((int)points[3].x, (int)points[3].y, (int)points[0].x, (int)points[0].y, color);
}

// グリッドの描画
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {

	const float kGridHalfWidth = 2.0f;                                      // Gridの半分の幅
	const uint32_t kSubdivision = 10;                                       // 分割数
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision); // 1つ分の長さ

	// 奥から手前への線を順々にに引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {

		float posX = -kGridHalfWidth + kGridEvery * xIndex;

		Vector3 start = { posX, 0.0f, -kGridHalfWidth };
		Vector3 end = { posX, 0.0f, kGridHalfWidth };

		start = Transform(start, Multiply(viewProjectionMatrix, viewportMatrix));
		end = Transform(end, Multiply(viewProjectionMatrix, viewportMatrix));

		// 左から右も同じように順々に引いていく
		for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
			float posZ = -kGridHalfWidth + kGridEvery * zIndex;

			Vector3 startZ = { -kGridHalfWidth, 0.0f, posZ };
			Vector3 endZ = { kGridHalfWidth, 0.0f, posZ };

			startZ = Transform(startZ, Multiply(viewProjectionMatrix, viewportMatrix));
			endZ = Transform(endZ, Multiply(viewProjectionMatrix, viewportMatrix));

			Novice::DrawLine((int)start.x, (int)start.y, (int)end.x, (int)end.y, 0x060606FF);
			Novice::DrawLine((int)startZ.x, (int)startZ.y, (int)endZ.x, (int)endZ.y, 0x060606FF);
		}
	}
}

// 球と平面の当たり判定
bool IsCollision(const Sphere& sphere, const Plane& plane) {

	// 平面と点の距離を求める
	float distance = Dot(plane.normal, sphere.center) - plane.distance;

	if (distance < 0) { distance *= -1.0f; }

	if (distance <= sphere.radius) { return true; }
	else { return false; }
}

// 正射影ベクトル
Vector3 Project(const Vector3& v1, const Vector3& v2) {

	return Multiply(Dot(v1, v2) / powf(Length(v2), 2), v2);
}

// 反射ベクトルを求める関数
Vector3 Reflect(const Vector3& input, const Vector3& normal) {

	Vector3 resultReflect;

	resultReflect = input - 2 * Dot(input, normal) * normal;

	return resultReflect;
}

const char kWindowTitle[] = "LE2B_01_アキモト_ワタル";
static const int kWindowWidth = 1280;
static const int kWindowHeight = 720;

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };
	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };

	Sphere sphere = {
		.center{ 0.0f, 0.0f, 0.0f },
		.radius{ 0.1f },
	};

	Plane plane = {
		.normal{ Normalize({ -0.2f, 0.9f, -0.3f })},
		.distance{ 0.0f },
	};

	Ball ball = {
		.position{ 0.8f, 1.2f, -0.3f },
		.acceleration{ 0.0f, -9.8f, 0.0f },
		.mass{ 2.0f },
		.radius{ 0.05f },
		.color{ WHITE },
	};

	float deltaTime = 1.0f / 60.0f;

	float e = 0.5f;

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

		ball.velocity = ball.velocity + ball.acceleration * deltaTime;
		ball.position = ball.position + ball.velocity * deltaTime;

		sphere.center = ball.position;

		if (IsCollision(sphere, plane)) {
			Vector3 reflected = Reflect(ball.velocity, plane.normal);
			Vector3 projectToNormal = Project(reflected, plane.normal);
			Vector3 movingDirection = reflected - projectToNormal;
			ball.velocity = projectToNormal * e + movingDirection;
		}

		ImGui::Begin("Window");
		ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::End();

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		DrawGrid(viewProjectionMatrix, viewportMatrix);

		DrawPlane(plane, viewProjectionMatrix, viewportMatrix, WHITE);

		DrawSphere(sphere, viewProjectionMatrix, viewportMatrix, WHITE);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
