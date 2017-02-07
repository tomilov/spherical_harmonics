//-----------------------------------------------------------------------------
// File: Framework\Math\Vector.h
// Copyright (c) 2007 Advanced Micro Devices Inc. All rights reserved.
//-----------------------------------------------------------------------------





#ifndef _VECTOR_H_
#define _VECTOR_H_

/** \file
	Standard 3D math definitions and routines
*/

#include "../Platform.h"
#include <math.h>

/** 2D single-precision floating-point vector */
struct float2
{
	float x, y;

	float2(){}
	float2(const float iv)
	{
		x = y = iv;
	}
	float2(const float ix, const float iy)
	{
		x = ix;
		y = iy;
	}

	operator float *(){ return &x; }
	operator const float *() const { return &x; }

	void operator += (const float2 &v);
	void operator -= (const float2 &v);
	void operator *= (const float s);
	void operator *= (const float2 &v);
	void operator /= (const float s);
	void operator /= (const float2 &v);
};

float2 operator + (const float2 &u, const float2 &v);
float2 operator + (const float2 &v, const float s);
float2 operator - (const float2 &u, const float2 &v);
float2 operator - (const float2 &v, const float s);
float2 operator - (const float2 &v);
float2 operator * (const float2 &u, const float2 &v);
float2 operator * (const float s, const float2 &v);
float2 operator * (const float2 &v, const float s);
float2 operator / (const float2 &u, const float2 &v);
float2 operator / (const float2 &v, const float s);

bool operator == (const float2 &u, const float2 &v);
bool operator != (const float2 &u, const float2 &v);

/** 3D single-precision floating-point vector */
struct float3
{
	float x, y, z;

	float3(){}
	float3(const float iv)
	{
		x = y = z = iv;
	}
	float3(const float2 &xy, const float iz)
	{
		x = xy.x;
		y = xy.y;
		z = iz;
	}
	float3(const float ix, const float iy, const float iz)
	{
		x = ix;
		y = iy;
		z = iz;
	}

	operator float *(){ return &x; }
	operator const float *() const { return &x; }

	float2 xy() const { return float2(x, y); }

	float3 xzy() const { return float3(x, z, y); }
	float3 zyx() const { return float3(z, y, x); }
	float3 yxz() const { return float3(y, x, z); }
	float3 yzx() const { return float3(y, z, x); }
	float3 zxy() const { return float3(z, x, y); }

	void operator += (const float3 &v);
	void operator -= (const float3 &v);
	void operator *= (const float s);
	void operator *= (const float3 &v);
	void operator /= (const float s);
	void operator /= (const float3 &v);
};

float3 operator + (const float3 &u, const float3 &v);
float3 operator + (const float3 &v, const float s);
float3 operator - (const float3 &u, const float3 &v);
float3 operator - (const float3 &v, const float s);
float3 operator - (const float3 &v);
float3 operator * (const float3 &u, const float3 &v);
float3 operator * (const float s, const float3 &v);
float3 operator * (const float3 &v, const float s);
float3 operator / (const float3 &u, const float3 &v);
float3 operator / (const float3 &v, const float s);

bool operator == (const float3 &u, const float3 &v);
bool operator != (const float3 &u, const float3 &v);

/** 4D single-precision floating-point vector */
struct float4
{
	float x, y, z, w;

	float4(){}
	float4(const float iv)
	{
		x = y = z = w = iv;
	}
	float4(const float ix, const float iy, const float iz, const float iw)
	{
		x = ix;
		y = iy;
		z = iz;
		w = iw;
	}
	float4(const float2 &xy, const float iz, const float iw)
	{
		x = xy.x;
		y = xy.y;
		z = iz;
		w = iw;
	}
	float4(const float2 &xy, const float2 &zw)
	{
		x = xy.x;
		y = xy.y;
		z = zw.x;
		w = zw.y;
	}
	float4(const float3 &xyz, const float iw)
	{
		x = xyz.x;
		y = xyz.y;
		z = xyz.z;
		w = iw;
	}

	operator float *(){ return &x; }
	operator const float *() const { return &x; }

	float2 xy()  const { return float2(x, y); }
	float3 xyz() const { return float3(x, y, z); }

	float3 xzy() const { return float3(x, z, y); }
	float3 zyx() const { return float3(z, y, x); }
	float3 yxz() const { return float3(y, x, z); }
	float3 yzx() const { return float3(y, z, x); }
	float3 zxy() const { return float3(z, x, y); }

	void operator += (const float4 &v);
	void operator -= (const float4 &v);
	void operator *= (const float s);
	void operator *= (const float4 &v);
	void operator /= (const float s);
	void operator /= (const float4 &v);
};

float4 operator + (const float4 &u, const float4 &v);
float4 operator + (const float4 &v, const float s);
float4 operator - (const float4 &u, const float4 &v);
float4 operator - (const float4 &v, const float s);
float4 operator - (const float4 &v);
float4 operator * (const float4 &u, const float4 &v);
float4 operator * (const float s, const float4 &v);
float4 operator * (const float4 &v, const float s);
float4 operator / (const float4 &u, const float4 &v);
float4 operator / (const float4 &v, const float s);

bool operator == (const float4 &u, const float4 &v);
bool operator != (const float4 &u, const float4 &v);


/** 2D half-precision floating-point vector */
struct half2
{
	half x, y;

	half2(){}
	half2(const half iv)
	{
		x = y = iv;
	}
	half2(const half ix, const half iy)
	{
		x = ix;
		y = iy;
	}
	half2(const float2 &iv)
	{
		x = iv.x;
		y = iv.y;
	}
	operator const half *() const { return &x; }
};

/** 3D half-precision floating-point vector */
struct half3
{
	half x, y, z;

	half3(){}
	half3(const half iv)
	{
		x = y = z = iv;
	}
	half3(const half ix, const half iy, const half iz)
	{
		x = ix;
		y = iy;
		z = iz;
	}
	half3(const float3 &iv)
	{
		x = iv.x;
		y = iv.y;
		z = iv.z;
	}
	operator const half *() const { return &x; }
};

/** 4D half-precision floating-point vector */
struct half4
{
	half x, y, z, w;

	half4(){}
	half4(const half iv)
	{
		x = y = z = w = iv;
	}
	half4(const half ix, const half iy, const half iz, const half iw)
	{
		x = ix;
		y = iy;
		z = iz;
		w = iw;
	}
	half4(const float4 &iv)
	{
		x = iv.x;
		y = iv.y;
		z = iv.z;
		w = iv.w;
	}
	operator const half *() const { return &x; }
};

/** 4D integer vector */
struct int4 {
	int x, y, z, w;

	int4(){}
	int4(const int iv)
	{
		x = y = z = w = iv;
	}
	int4(const int ix, const int iy, const int iz, const int iw)
	{
		x = ix;
		y = iy;
		z = iz;
		w = iw;
	}
	operator int *(){ return &x; }
	operator const int *() const { return &x; }
};


float dot(const float2 &u, const float2 &v);
float dot(const float3 &u, const float3 &v);
float dot(const float4 &u, const float4 &v);

float lerp(const float u, const float v, const float x);
float2 lerp(const float2 &u, const float2 &v, const float x);
float3 lerp(const float3 &u, const float3 &v, const float x);
float4 lerp(const float4 &u, const float4 &v, const float x);

float cerp(const float u0, const float u1, const float u2, const float u3, const float x);
float2 cerp(const float2 &u0, const float2 &u1, const float2 &u2, const float2 &u3, const float x);
float3 cerp(const float3 &u0, const float3 &u1, const float3 &u2, const float3 &u3, const float x);
float4 cerp(const float4 &u0, const float4 &u1, const float4 &u2, const float4 &u3, const float x);

float herp(const float u0, const float u1, const float u2, const float u3, const float x, const float tension, const float bias);
float2 herp(const float2 u0, const float2 u1, const float2 u2, const float2 u3, const float x, const float tension, const float bias);
float3 herp(const float3 u0, const float3 u1, const float3 u2, const float3 u3, const float x, const float tension, const float bias);
float4 herp(const float4 u0, const float4 u1, const float4 u2, const float4 u3, const float x, const float tension, const float bias);

/** Fit t to an S-Curve */
inline float sCurve(const float t)
{
	return t * t * (3 - 2 * t);
}


#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

/** Returns the lesser of two vectors */
template <typename DATA_TYPE>
inline DATA_TYPE min(const DATA_TYPE x, const DATA_TYPE y)
{
	return (x < y)? x : y;
}

/** Returns the greater of two vectors */
template <typename DATA_TYPE>
inline DATA_TYPE max(const DATA_TYPE x, const DATA_TYPE y)
{
	return (x > y)? x : y;
}

float2 min(const float2 u, const float2 v);
float2 max(const float2 u, const float2 v);
float3 min(const float3 u, const float3 v);
float3 max(const float3 u, const float3 v);
float4 min(const float4 u, const float4 v);
float4 max(const float4 u, const float4 v);

/** Clamps x to [lower-upper] */
template <typename DATA_TYPE>
inline DATA_TYPE clamp(const DATA_TYPE x, const float lower, const float upper)
{
	return max(min(x, DATA_TYPE(upper)), DATA_TYPE(lower));
}

/** Clamps x to [0-1] */
#define saturate(x) clamp(x, 0, 1)

int signP(const int x);     /**< Returns +1 if x > 0, -1 otherwise */
int signN(const int x);     /**< Returns +1 if x < 0, -1 otherwise */
int sign(const int x);      
float sign(const float x);
float2 sign(const float2 &v);
float3 sign(const float3 &v);
float4 sign(const float4 &v);

float smoothstep(const float edge0, const float edge1, const float x);

/** Compute the square of the input x */
template <typename DATA_TYPE>
inline DATA_TYPE sqr(const DATA_TYPE x)
{
	return x * x;
}

float2 normalize(const float2 &v);
float3 normalize(const float3 &v);
float4 normalize(const float4 &v);

float length(const float2 &v);
float length(const float3 &v);
float length(const float4 &v);

/** Compute the distance between two nD vectors, where n is [2-4] */
//#define distance(u, v) length((u) - (v))
template< typename T >
float distance(T u, T v)
{
    return length(u - v);
}

/** Compute the cross product of two 3D vectors */
float3 cross(const float3 &u, const float3 &v);

/** Convert RGBE format to RGB format */
float3 rgbeToRGB(unsigned char *rgbe);

/** Convert 4 floats into a single DWORD value. @{ */
unsigned int toRGBA(const float4 &u);
unsigned int toBGRA(const float4 &u);
//@}





/** 2x2 matrix */
struct float2x2
{
	float elem[2][2];

	float2x2(){}
	float2x2(const float m00, const float m01,
	         const float m10, const float m11)
	{
		elem[0][0] = m00; elem[0][1] = m01;
		elem[1][0] = m10; elem[1][1] = m11;
	}
	float2x2(const float2 &v0, const float2 &v1)
	{
		elem[0][0] = v0.x; elem[0][1] = v0.y;
		elem[1][0] = v1.x; elem[1][1] = v1.y;
	}

	operator const float *() const { return (const float *) elem; }
};

float2x2 operator + (const float2x2 &m, const float2x2 &n);
float2x2 operator - (const float2x2 &m, const float2x2 &n);
float2x2 operator - (const float2x2 &m);
float2x2 operator * (const float2x2 &m, const float2x2 &n);
float2 operator * (const float2x2 &m, const float2 &v);
float2x2 operator * (const float2x2 &m, const float x);

float2x2 transpose(const float2x2 &m);

/** 3x3 matrix */
struct float3x3
{
	float elem[3][3];

	float3x3(){}
	float3x3(const float m00, const float m01, const float m02,
	         const float m10, const float m11, const float m12,
	         const float m20, const float m21, const float m22)
	{
		elem[0][0] = m00; elem[0][1] = m01; elem[0][2] = m02;
		elem[1][0] = m10; elem[1][1] = m11; elem[1][2] = m12;
		elem[2][0] = m20; elem[2][1] = m21; elem[2][2] = m22;
	}
	float3x3(const float3 &v0, const float3 &v1, const float3 &v2)
	{
		elem[0][0] = v0.x; elem[0][1] = v0.y; elem[0][2] = v0.z;
		elem[1][0] = v1.x; elem[1][1] = v1.y; elem[1][2] = v1.z;
		elem[2][0] = v2.x; elem[2][1] = v2.y; elem[2][2] = v2.z;
	}

	float3 &getRow(const unsigned int index){ return *((float3 *) elem[index]); }

	operator const float *() const { return (const float *) elem; }
};

float3x3 operator + (const float3x3 &m, const float3x3 &n);
float3x3 operator - (const float3x3 &m, const float3x3 &n);
float3x3 operator - (const float3x3 &m);
float3x3 operator * (const float3x3 &m, const float3x3 &n);
float3 operator * (const float3x3 &m, const float3 &v);
float3x3 operator * (const float3x3 &m, const float x);

float3x3 transpose(const float3x3 &m);
float3x3 operator ! (const float3x3 &m);

/** 4x4 matrix */
struct float4x4
{
	float elem[4][4];

	float4x4(){}
	float4x4(const float m00, const float m01, const float m02, const float m03,
	         const float m10, const float m11, const float m12, const float m13,
	         const float m20, const float m21, const float m22, const float m23,
	         const float m30, const float m31, const float m32, const float m33)
	{
		elem[0][0] = m00; elem[0][1] = m01; elem[0][2] = m02; elem[0][3] = m03;
		elem[1][0] = m10; elem[1][1] = m11; elem[1][2] = m12; elem[1][3] = m13;
		elem[2][0] = m20; elem[2][1] = m21; elem[2][2] = m22; elem[2][3] = m23;
		elem[3][0] = m30; elem[3][1] = m31; elem[3][2] = m32; elem[3][3] = m33;
	}
	float4x4(const float4 &v0, const float4 &v1, const float4 &v2, const float4 &v3)
	{
		elem[0][0] = v0.x; elem[0][1] = v0.y; elem[0][2] = v0.z; elem[0][3] = v0.w;
		elem[1][0] = v1.x; elem[1][1] = v1.y; elem[1][2] = v1.z; elem[1][3] = v1.w;
		elem[2][0] = v2.x; elem[2][1] = v2.y; elem[2][2] = v2.z; elem[2][3] = v2.w;
		elem[3][0] = v3.x; elem[3][1] = v3.y; elem[3][2] = v3.z; elem[3][3] = v3.w;
	}

	float3 getRightVec()   const { return float3(elem[0][0], elem[0][1], elem[0][2]); }
	float3 getUpVec()      const { return float3(elem[1][0], elem[1][1], elem[1][2]); }
	float3 getForwardVec() const { return float3(elem[2][0], elem[2][1], elem[2][2]); }

	operator const float *() const { return (const float *) elem; }
	void translate(const float3 &v);
};

float4x4 operator + (const float4x4 &m, const float4x4 &n);
float4x4 operator - (const float4x4 &m, const float4x4 &n);
float4x4 operator - (const float4x4 &m);
float4x4 operator * (const float4x4 &m, const float4x4 &n);
float4 operator * (const float4x4 &m, const float4 &v);
float4x4 operator * (const float4x4 &m, const float x);

float4x4 transpose(const float4x4 &m);
float4x4 operator ! (const float4x4 &m);

/** Computes a 2x2 rotation matrix about a given angle */
float2x2 rotate2(const float angle);

/** Computes an NxN rotation matrix about a (set of) given angle(s). @{ */
float3x3 rotateX3(const float angle);
float3x3 rotateY3(const float angle);
float3x3 rotateZ3(const float angle);
float3x3 rotateXY3(const float angleX, const float angleY);
float3x3 rotateZXY3(const float angleX, const float angleY, const float angleZ);

float4x4 rotateX4(const float angle);
float4x4 rotateY4(const float angle);
float4x4 rotateZ4(const float angle);
float4x4 rotateXY4(const float angleX, const float angleY);
float4x4 rotateZXY4(const float angleX, const float angleY, const float angleZ);
//@}

/** Returns an NxN scaling matrix. @{ */
float2x2 scale2(const float x, const float y);
float3x3 scale3(const float x, const float y, const float z);
float4x4 scale4(const float x, const float y, const float z);
//@}
/** Creates a 4x4 translation matrix by the input vector */
float4x4 translate(const float x, const float y, const float z);
float4x4 translate(const float3 &v);

/** Creates an orthographic projection matrix */
float4x4 OrthoMatrix(const float left, const float right, const float top, const float bottom, const float zNear, const float zFar);

/** Creates a perspective-correct projection matrix */
float4x4 PerspectiveMatrixX(const float fov, const int width, const int height, const float zNear, const float zFar);
float4x4 PerspectiveMatrixY(const float fov, const int width, const int height, const float zNear, const float zFar);

/** Creates a shadow projection matrix */
float4x4 ShadowMatrix(const float3 &planeNormal, const float planeOffset, const float3 &lightPos);

/** Creates a mirror matrix for a given plane */
float4x4 MirrorMatrix(const float3 &planeNormal, const float planeOffset);

/** Creates a modelview matrix for rendering to a given cubemap face */
float4x4 CubemapModelviewMatrix(const unsigned int face);

/** Creates a projection matrix for rendering to a cubemap */
float4x4 CubemapProjectionMatrix(const float zNear, const float zFar);

/** Creates a standard view matrix */
float4x4 LookAt(const float3 &from, const float3 &at);

/** Converts the z projection from GL to D3D convention */
float4x4 OpenGLToD3DProjectionMatrix(const float4x4 &m);

/** Pegs the geometry to the far clipping plane. Useful for instance for skyboxes. */
float4x4 PegToFarPlane(const float4x4 &m);

/** Returns an NxN identity matrix. @{ */
float2x2 Identity2();
float3x3 Identity3();
float4x4 Identity4();
//@}

#endif
