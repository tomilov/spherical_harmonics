//-----------------------------------------------------------------------------
// File: Framework\Math\Vector.cpp
// Copyright (c) 2007 Advanced Micro Devices Inc. All rights reserved.
//-----------------------------------------------------------------------------





#include "Vector.h"
#include <memory.h>

void float2::operator += (const float2 &v)
{
	x += v.x;
	y += v.y;
}

void float2::operator -= (const float2 &v)
{
	x -= v.x;
	y -= v.y;
}

void float2::operator *= (const float s)
{
	x *= s;
	y *= s;
}

void float2::operator *= (const float2 &v)
{
	x *= v.x;
	y *= v.y;
}

void float2::operator /= (const float s)
{
	x /= s;
	y /= s;
}

void float2::operator /= (const float2 &v)
{
	x /= v.x;
	y /= v.y;
}

float2 operator + (const float2 &u, const float2 &v)
{
	return float2(u.x + v.x, u.y + v.y);
}

float2 operator + (const float2 &v, const float s)
{
	return float2(v.x + s, v.y + s);
}

float2 operator - (const float2 &u, const float2 &v)
{
	return float2(u.x - v.x, u.y - v.y);
}

float2 operator - (const float2 &v, const float s)
{
	return float2(v.x - s, v.y - s);
}

float2 operator - (const float2 &v)
{
	return float2(-v.x, -v.y);
}

float2 operator * (const float2 &u, const float2 &v)
{
	return float2(u.x * v.x, u.y * v.y);
}

float2 operator * (const float s, const float2 &v)
{
	return float2(v.x * s, v.y * s);
}

float2 operator * (const float2 &v, const float s)
{
	return float2(v.x * s, v.y * s);
}

float2 operator / (const float2 &u, const float2 &v)
{
	return float2(u.x / v.x, u.y / v.y);
}

float2 operator / (const float2 &v, const float s)
{
	return float2(v.x / s, v.y / s);
}

bool operator == (const float2 &u, const float2 &v)
{
    return (u.x == v.x && u.y == v.y);
}

bool operator != (const float2 &u, const float2 &v)
{
    return (u.x != v.x || u.y != v.y);
}


void float3::operator += (const float3 &v)
{
	x += v.x;
	y += v.y;
	z += v.z;
}

void float3::operator -= (const float3 &v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
}

void float3::operator *= (const float s)
{
	x *= s;
	y *= s;
	z *= s;
}

void float3::operator *= (const float3 &v)
{
	x *= v.x;
	y *= v.y;
	z *= v.z;
}

void float3::operator /= (const float s)
{
	x /= s;
	y /= s;
	z /= s;
}

void float3::operator /= (const float3 &v)
{
	x /= v.x;
	y /= v.y;
	z /= v.z;
}

float3 operator + (const float3 &u, const float3 &v)
{
	return float3(u.x + v.x, u.y + v.y, u.z + v.z);
}

float3 operator + (const float3 &v, const float s)
{
	return float3(v.x + s, v.y + s, v.z + s);
}

float3 operator - (const float3 &u, const float3 &v)
{
	return float3(u.x - v.x, u.y - v.y, u.z - v.z);
}

float3 operator - (const float3 &v, const float s)
{
	return float3(v.x - s, v.y - s, v.z - s);
}

float3 operator - (const float3 &v)
{
	return float3(-v.x, -v.y, -v.z);
}

float3 operator * (const float3 &u, const float3 &v)
{
	return float3(u.x * v.x, u.y * v.y, u.z * v.z);
}

float3 operator * (const float s, const float3 &v)
{
	return float3(v.x * s, v.y * s, v.z * s);
}

float3 operator * (const float3 &v, const float s)
{
	return float3(v.x * s, v.y * s, v.z * s);
}

float3 operator / (const float3 &u, const float3 &v)
{
	return float3(u.x / v.x, u.y / v.y, u.z / v.z);
}

float3 operator / (const float3 &v, const float s)
{
	return float3(v.x / s, v.y / s, v.z / s);
}

bool operator == (const float3 &u, const float3 &v)
{
    return (u.x == v.x && u.y == v.y && u.z == v.z);
}

bool operator != (const float3 &u, const float3 &v)
{
    return (u.x != v.x || u.y != v.y || u.z != v.z);
}


void float4::operator += (const float4 &v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	w += v.w;
}

void float4::operator -= (const float4 &v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	w -= v.w;
}

void float4::operator *= (const float s)
{
	x *= s;
	y *= s;
	z *= s;
	w *= s;
}

void float4::operator *= (const float4 &v)
{
	x *= v.x;
	y *= v.y;
	z *= v.z;
	w *= v.w;
}

void float4::operator /= (const float s)
{
	x /= s;
	y /= s;
	z /= s;
	w /= s;
}

void float4::operator /= (const float4 &v)
{
	x /= v.x;
	y /= v.y;
	z /= v.z;
	w /= v.w;
}

float4 operator + (const float4 &u, const float4 &v)
{
	return float4(u.x + v.x, u.y + v.y, u.z + v.z, u.w + v.w);
}

float4 operator + (const float4 &v, const float s)
{
	return float4(v.x + s, v.y + s, v.z + s, v.w + s);
}

float4 operator - (const float4 &u, const float4 &v)
{
	return float4(u.x - v.x, u.y - v.y, u.z - v.z, u.w - v.w);
}

float4 operator - (const float4 &v, const float s)
{
	return float4(v.x - s, v.y - s, v.z - s, v.w - s);
}

float4 operator - (const float4 &v)
{
	return float4(-v.x, -v.y, -v.z, -v.w);
}

float4 operator * (const float4 &u, const float4 &v)
{
	return float4(u.x * v.x, u.y * v.y, u.z * v.z, u.w * v.w);
}

float4 operator * (const float s, const float4 &v)
{
	return float4(v.x * s, v.y * s, v.z * s, v.w * s);
}

float4 operator * (const float4 &v, const float s)
{
	return float4(v.x * s, v.y * s, v.z * s, v.w * s);
}

float4 operator / (const float4 &u, const float4 &v)
{
	return float4(u.x / v.x, u.y / v.y, u.z / v.z, u.w / v.w);
}

float4 operator / (const float4 &v, const float s)
{
	return float4(v.x / s, v.y / s, v.z / s, v.w / s);
}

bool operator == (const float4 &u, const float4 &v)
{
    return (u.x == v.x && u.y == v.y && u.z == v.z && u.w == v.w);
}

bool operator != (const float4 &u, const float4 &v)
{
    return (u.x != v.x || u.y != v.y || u.z != v.z || u.w != v.w);
}




float dot(const float2 &u, const float2 &v)
{
	return u.x * v.x + u.y * v.y;
}

float dot(const float3 &u, const float3 &v
		  ){
	return u.x * v.x + u.y * v.y + u.z * v.z;
}

float dot(const float4 &u, const float4 &v)
{
	return u.x * v.x + u.y * v.y + u.z * v.z + u.w * v.w;
}

float lerp(const float u, const float v, const float x)
{
	return u * (1 - x) + v * x;
}

float2 lerp(const float2 &u, const float2 &v, const float x)
{
	return u * (1 - x) + v * x;
}

float3 lerp(const float3 &u, const float3 &v, const float x)
{
	return u * (1 - x) + v * x;
}

float4 lerp(const float4 &u, const float4 &v, const float x)
{
	return u * (1 - x) + v * x;
}

float cerp(const float u0, const float u1, const float u2, const float u3, const float x)
{
	float p = (u3 - u2) - (u0 - u1);
	float q = (u0 - u1) - p;
	float r = u2 - u0;
	return x * (x * (x * p + q) + r) + u1;
}

float2 cerp(const float2 &u0, const float2 &u1, const float2 &u2, const float2 &u3, const float x)
{
	float2 p = (u3 - u2) - (u0 - u1);
	float2 q = (u0 - u1) - p;
	float2 r = u2 - u0;
	return x * (x * (x * p + q) + r) + u1;
}

float3 cerp(const float3 &u0, const float3 &u1, const float3 &u2, const float3 &u3, const float x)
{
	float3 p = (u3 - u2) - (u0 - u1);
	float3 q = (u0 - u1) - p;
	float3 r = u2 - u0;
	return x * (x * (x * p + q) + r) + u1;
}

float4 cerp(const float4 &u0, const float4 &u1, const float4 &u2, const float4 &u3, const float x)
{
	float4 p = (u3 - u2) - (u0 - u1);
	float4 q = (u0 - u1) - p;
	float4 r = u2 - u0;
	return x * (x * (x * p + q) + r) + u1;
}

float herp(const float u0, const float u1, const float u2, const float u3, const float x, const float tension, const float bias)
{
	float m0 = (u1 - u0) * (1 + bias) * (1 - tension) / 2 + (u2 - u1) * (1 - bias) * (1 - tension) / 2;
	float m1 = (u2 - u1) * (1 + bias) * (1 - tension) / 2 + (u3 - u2) * (1 - bias) * (1 - tension) / 2;
	float x2 = x * x;
	float x3 = x * x2;
	float a0 =  2 * x3 - 3 * x2 + 1;
	float a1 =      x3 - 2 * x2 + x;
	float a2 =      x3 -     x2;
	float a3 = -2 * x3 + 3 * x2;

	return a0 * u1 + a1 * m0 + a2 * m1 + a3 * u2;
}

float2 herp(const float2 u0, const float2 u1, const float2 u2, const float2 u3, const float x, const float tension, const float bias)
{
	float2 m0 = (u1 - u0) * (1 + bias) * (1 - tension) / 2 + (u2 - u1) * (1 - bias) * (1 - tension) / 2;
	float2 m1 = (u2 - u1) * (1 + bias) * (1 - tension) / 2 + (u3 - u2) * (1 - bias) * (1 - tension) / 2;
	float x2 = x * x;
	float x3 = x * x2;
	float a0 =  2 * x3 - 3 * x2 + 1;
	float a1 =      x3 - 2 * x2 + x;
	float a2 =      x3 -     x2;
	float a3 = -2 * x3 + 3 * x2;

	return a0 * u1 + a1 * m0 + a2 * m1 + a3 * u2;
}

float3 herp(const float3 u0, const float3 u1, const float3 u2, const float3 u3, const float x, const float tension, const float bias)
{
	float3 m0 = (u1 - u0) * (1 + bias) * (1 - tension) / 2 + (u2 - u1) * (1 - bias) * (1 - tension) / 2;
	float3 m1 = (u2 - u1) * (1 + bias) * (1 - tension) / 2 + (u3 - u2) * (1 - bias) * (1 - tension) / 2;
	float x2 = x * x;
	float x3 = x * x2;
	float a0 =  2 * x3 - 3 * x2 + 1;
	float a1 =      x3 - 2 * x2 + x;
	float a2 =      x3 -     x2;
	float a3 = -2 * x3 + 3 * x2;

	return a0 * u1 + a1 * m0 + a2 * m1 + a3 * u2;
}

float4 herp(const float4 u0, const float4 u1, const float4 u2, const float4 u3, const float x, const float tension, const float bias)
{
	float4 m0 = (u1 - u0) * (1 + bias) * (1 - tension) / 2 + (u2 - u1) * (1 - bias) * (1 - tension) / 2;
	float4 m1 = (u2 - u1) * (1 + bias) * (1 - tension) / 2 + (u3 - u2) * (1 - bias) * (1 - tension) / 2;
	float x2 = x * x;
	float x3 = x * x2;
	float a0 =  2 * x3 - 3 * x2 + 1;
	float a1 =      x3 - 2 * x2 + x;
	float a2 =      x3 -     x2;
	float a3 = -2 * x3 + 3 * x2;

	return a0 * u1 + a1 * m0 + a2 * m1 + a3 * u2;
}



float2 min(const float2 u, const float2 v)
{
    return float2((u.x < v.x)? u.x : v.x, (u.y < v.y)? u.y : v.y);
}

float2 max(const float2 u, const float2 v)
{
    return float2((u.x > v.x)? u.x : v.x, (u.y > v.y)? u.y : v.y);
}

float3 min(const float3 u, const float3 v)
{
    return float3((u.x < v.x)? u.x : v.x, (u.y < v.y)? u.y : v.y, (u.z < v.z)? u.z : v.z);
}

float3 max(const float3 u, const float3 v)
{
    return float3((u.x > v.x)? u.x : v.x, (u.y > v.y)? u.y : v.y, (u.z > v.z)? u.z : v.z);
}

float4 min(const float4 u, const float4 v)
{
    return float4((u.x < v.x)? u.x : v.x, (u.y < v.y)? u.y : v.y, (u.z < v.z)? u.z : v.z, (u.w < v.w)? u.w : v.w);
}

float4 max(const float4 u, const float4 v)
{
    return float4((u.x > v.x)? u.x : v.x, (u.y > v.y)? u.y : v.y, (u.z > v.z)? u.z : v.z, (u.w > v.w)? u.w : v.w);
}

int signP(const int x)
{
    return 1 - 2 * (x < 0);
}

int signN(const int x)
{
    return 2 * (x > 0) - 1;
}

int sign(const int x)
{
    return (x > 0) - (x < 0);
}

float sign(const float x)
{
	return (x < 0)? -1.0f : 1.0f;
}

float2 sign(const float2 &v)
{
	return float2((v.x < 0)? -1.0f : 1.0f, (v.y < 0)? -1.0f : 1.0f);
}

float3 sign(const float3 &v)
{
	return float3((v.x < 0)? -1.0f : 1.0f, (v.y < 0)? -1.0f : 1.0f, (v.z < 0)? -1.0f : 1.0f);
}

float4 sign(const float4 &v)
{
	return float4((v.x < 0)? -1.0f : 1.0f, (v.y < 0)? -1.0f : 1.0f, (v.z < 0)? -1.0f : 1.0f, (v.w < 0)? -1.0f : 1.0f);
}

float smoothstep(const float edge0, const float edge1, const float x)
{
	float t = saturate((x - edge0) / (edge1 - edge0));
	return t * t * (3 - 2 * t);
}

float2 normalize(const float2 &v)
{
	float invLen = 1.0f / sqrtf(v.x * v.x + v.y * v.y);
	return v * invLen;
}

float3 normalize(const float3 &v)
{
	float invLen = 1.0f / sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	return v * invLen;
}

float4 normalize(const float4 &v)
{
	float invLen = 1.0f / sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
	return v * invLen;
}

float length(const float2 &v)
{
	return sqrtf(v.x * v.x + v.y * v.y);
}

float length(const float3 &v)
{
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

float length(const float4 &v)
{
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
}

float3 cross(const float3 &u, const float3 &v)
{
	return float3(u.y * v.z - v.y * u.z, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
}

float3 rgbeToRGB(unsigned char *rgbe)
{
	if (rgbe[3]){
		return float3(rgbe[0], rgbe[1], rgbe[2]) * ldexpf(1.0f, rgbe[3] - (int) (128 + 8));
	} else return float3(0,0,0);
}

unsigned int toRGBA(const float4 &u)
{
	return (int(u.x * 255) | (int(u.y * 255) << 8) | (int(u.z * 255) << 16) | (int(u.w * 255) << 24));
}

unsigned int toBGRA(const float4 &u)
{
	return (int(u.z * 255) | (int(u.y * 255) << 8) | (int(u.x * 255) << 16) | (int(u.w * 255) << 24));
}




float2x2 operator + (const float2x2 &m, const float2x2 &n)
{
	return float2x2(
		m.elem[0][0] + n.elem[0][0], m.elem[0][1] + n.elem[0][1],
		m.elem[1][0] + n.elem[1][0], m.elem[1][1] + n.elem[1][1]);
}

float2x2 operator - (const float2x2 &m, const float2x2 &n)
{
	return float2x2(
		m.elem[0][0] - n.elem[0][0], m.elem[0][1] - n.elem[0][1],
		m.elem[1][0] - n.elem[1][0], m.elem[1][1] - n.elem[1][1]);
}

float2x2 operator - (const float2x2 &m)
{
	return float2x2(-m.elem[0][0], -m.elem[0][1], -m.elem[1][0], -m.elem[1][1]);
}

float2x2 operator * (const float2x2 &m, const float2x2 &n)
{
	return float2x2(
		m.elem[0][0] * n.elem[0][0] + m.elem[0][1] * n.elem[1][0],   m.elem[0][0] * n.elem[0][1] + m.elem[0][1] * n.elem[1][1],
		m.elem[1][0] * n.elem[0][0] + m.elem[1][1] * n.elem[1][0],   m.elem[1][0] * n.elem[0][1] + m.elem[1][1] * n.elem[1][1]);
}

float2 operator * (const float2x2 &m, const float2 &v)
{
	return float2(m.elem[0][0] * v.x + m.elem[0][1] * v.y, m.elem[1][0] * v.x + m.elem[1][1] * v.y);
}

float2x2 operator * (const float2x2 &m, const float x)
{
	float2x2 mat;
	mat.elem[0][0] = m.elem[0][0] * x;
	mat.elem[0][1] = m.elem[0][1] * x;
	mat.elem[1][0] = m.elem[1][0] * x;
	mat.elem[1][1] = m.elem[1][1] * x;
	return mat;
}

float2x2 transpose(const float2x2 &m)
{
	return float2x2(
		m.elem[0][0], m.elem[1][0],
		m.elem[0][1], m.elem[1][1]);
}


float3x3 operator + (const float3x3 &m, const float3x3 &n)
{
	return float3x3(
		m.elem[0][0] + n.elem[0][0], m.elem[0][1] + n.elem[0][1], m.elem[0][2] + n.elem[0][2],
		m.elem[1][0] + n.elem[1][0], m.elem[1][1] + n.elem[1][1], m.elem[1][2] + n.elem[1][2],
		m.elem[2][0] + n.elem[2][0], m.elem[2][1] + n.elem[2][1], m.elem[2][2] + n.elem[2][2]);
}

float3x3 operator - (const float3x3 &m, const float3x3 &n)
{
	return float3x3(
		m.elem[0][0] - n.elem[0][0], m.elem[0][1] - n.elem[0][1], m.elem[0][2] - n.elem[0][2],
		m.elem[1][0] - n.elem[1][0], m.elem[1][1] - n.elem[1][1], m.elem[1][2] - n.elem[1][2],
		m.elem[2][0] - n.elem[2][0], m.elem[2][1] - n.elem[2][1], m.elem[2][2] - n.elem[2][2]);
}

float3x3 operator - (const float3x3 &m)
{
	return float3x3(
		-m.elem[0][0], -m.elem[0][1], -m.elem[0][2],
		-m.elem[1][0], -m.elem[1][1], -m.elem[1][2],
		-m.elem[2][0], -m.elem[2][1], -m.elem[2][2]);
}

float3x3 operator * (const float3x3 &m, const float3x3 &n)
{
	float3x3 mat;

	for (unsigned int k = 0; k < 3; k++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			mat.elem[k][j] = 0;
			for (unsigned int i = 0; i < 3; i++)
			{
				mat.elem[k][j] += m.elem[k][i] * n.elem[i][j];
			}
		}
	}
	return mat;
}

float3 operator * (const float3x3 &m, const float3 &v)
{
	return float3(
		m.elem[0][0] * v.x + m.elem[0][1] * v.y + m.elem[0][2] * v.z,
		m.elem[1][0] * v.x + m.elem[1][1] * v.y + m.elem[1][2] * v.z,
		m.elem[2][0] * v.x + m.elem[2][1] * v.y + m.elem[2][2] * v.z);
}

float3x3 operator * (const float3x3 &m, const float x)
{
	float3x3 mat;
	for (unsigned int j = 0; j < 3; j++)
	{
		for (unsigned int i = 0; i < 3; i++)
		{
			mat.elem[j][i] = m.elem[j][i] * x;
		}
	}
	return mat;
}

float3x3 transpose(const float3x3 &m)
{
	return float3x3(
		m.elem[0][0], m.elem[1][0], m.elem[2][0],
		m.elem[0][1], m.elem[1][1], m.elem[2][1],
		m.elem[0][2], m.elem[1][2], m.elem[2][2]);
}

float3x3 operator ! (const float3x3 &m)
{
	float p0 = m.elem[1][1] * m.elem[2][2] - m.elem[2][1] * m.elem[1][2];
	float p1 = m.elem[1][2] * m.elem[2][0] - m.elem[1][0] * m.elem[2][2];
	float p2 = m.elem[1][0] * m.elem[2][1] - m.elem[1][1] * m.elem[2][0];

	float invDet = 1.0f / (m.elem[0][0] * p0 + m.elem[0][1] * p1 + m.elem[0][2] * p2);

	return float3x3(
		p0, m.elem[2][1] * m.elem[0][2] - m.elem[0][1] * m.elem[2][2], m.elem[0][1] * m.elem[1][2] - m.elem[1][1] * m.elem[0][2],
		p1, m.elem[0][0] * m.elem[2][2] - m.elem[2][0] * m.elem[0][2], m.elem[1][0] * m.elem[0][2] - m.elem[0][0] * m.elem[1][2],
		p2, m.elem[2][0] * m.elem[0][1] - m.elem[0][0] * m.elem[2][1], m.elem[0][0] * m.elem[1][1] - m.elem[0][1] * m.elem[1][0]) * invDet;
}


void float4x4::translate(const float3 &v)
{
	elem[0][3] += (v.x * elem[0][0] + v.y * elem[0][1] + v.z * elem[0][2]);
	elem[1][3] += (v.x * elem[1][0] + v.y * elem[1][1] + v.z * elem[1][2]);
	elem[2][3] += (v.x * elem[2][0] + v.y * elem[2][1] + v.z * elem[2][2]);
	elem[3][3] += (v.x * elem[3][0] + v.y * elem[3][1] + v.z * elem[3][2]);
}

float4x4 operator + (const float4x4 &m, const float4x4 &n)
{
	float4x4 mat;
	for (unsigned int i = 0; i < 16; i++)
	{
		((float *) &mat)[i] = ((float *) &m)[i] + ((float *) &n)[i];
	}
	return mat;
}

float4x4 operator - (const float4x4 &m, const float4x4 &n)
{
	float4x4 mat;
	for (unsigned int i = 0; i < 16; i++)
	{
		((float *) &mat)[i] = ((float *) &m)[i] - ((float *) &n)[i];
	}
	return mat;
}

float4x4 operator - (const float4x4 &m)
{
	float4x4 mat;
	for (unsigned int i = 0; i < 16; i++)
	{
		((float *) &mat)[i] = -((float *) &m)[i];
	}
	return mat;
}

float4x4 operator * (const float4x4 &m, const float4x4 &n)
{
	float4x4 mat;
	for (unsigned int k = 0; k < 4; k++)
	{
		for (unsigned int j = 0; j < 4; j++)
		{
			mat.elem[k][j] = 0;
			for (unsigned int i = 0; i < 4; i++)
			{
				mat.elem[k][j] += m.elem[k][i] * n.elem[i][j];
			}
		}
	}
	return mat;
}

float4 operator * (const float4x4 &m, const float4 &v)
{
	return float4(
		m.elem[0][0] * v.x + m.elem[0][1] * v.y + m.elem[0][2] * v.z + m.elem[0][3] * v.w,
		m.elem[1][0] * v.x + m.elem[1][1] * v.y + m.elem[1][2] * v.z + m.elem[1][3] * v.w,
		m.elem[2][0] * v.x + m.elem[2][1] * v.y + m.elem[2][2] * v.z + m.elem[2][3] * v.w,
		m.elem[3][0] * v.x + m.elem[3][1] * v.y + m.elem[3][2] * v.z + m.elem[3][3] * v.w);
}

float4x4 operator * (const float4x4 &m, const float x)
{
	float4x4 mat;
	for (unsigned int j = 0; j < 4; j++)
	{
		for (unsigned int i = 0; i < 4; i++)
		{
			mat.elem[j][i] = m.elem[j][i] * x;
		}
	}
	return mat;
}

float4x4 transpose(const float4x4 &m)
{
	return float4x4(
		m.elem[0][0], m.elem[1][0], m.elem[2][0], m.elem[3][0],
		m.elem[0][1], m.elem[1][1], m.elem[2][1], m.elem[3][1],
		m.elem[0][2], m.elem[1][2], m.elem[2][2], m.elem[3][2],
		m.elem[0][3], m.elem[1][3], m.elem[2][3], m.elem[3][3]);
}

float4x4 operator ! (const float4x4 &m)
{
	float4x4 mat;

	float p00 = m.elem[2][2] * m.elem[3][3];
	float p01 = m.elem[3][2] * m.elem[2][3];
	float p02 = m.elem[1][2] * m.elem[3][3];
	float p03 = m.elem[3][2] * m.elem[1][3];
	float p04 = m.elem[1][2] * m.elem[2][3];
	float p05 = m.elem[2][2] * m.elem[1][3];
	float p06 = m.elem[0][2] * m.elem[3][3];
	float p07 = m.elem[3][2] * m.elem[0][3];
	float p08 = m.elem[0][2] * m.elem[2][3];
	float p09 = m.elem[2][2] * m.elem[0][3];
	float p10 = m.elem[0][2] * m.elem[1][3];
	float p11 = m.elem[1][2] * m.elem[0][3];

	mat.elem[0][0] = (p00 * m.elem[1][1] + p03 * m.elem[2][1] + p04 * m.elem[3][1]) - (p01 * m.elem[1][1] + p02 * m.elem[2][1] + p05 * m.elem[3][1]);
	mat.elem[0][1] = (p01 * m.elem[0][1] + p06 * m.elem[2][1] + p09 * m.elem[3][1]) - (p00 * m.elem[0][1] + p07 * m.elem[2][1] + p08 * m.elem[3][1]);
	mat.elem[0][2] = (p02 * m.elem[0][1] + p07 * m.elem[1][1] + p10 * m.elem[3][1]) - (p03 * m.elem[0][1] + p06 * m.elem[1][1] + p11 * m.elem[3][1]);
	mat.elem[0][3] = (p05 * m.elem[0][1] + p08 * m.elem[1][1] + p11 * m.elem[2][1]) - (p04 * m.elem[0][1] + p09 * m.elem[1][1] + p10 * m.elem[2][1]);
	mat.elem[1][0] = (p01 * m.elem[1][0] + p02 * m.elem[2][0] + p05 * m.elem[3][0]) - (p00 * m.elem[1][0] + p03 * m.elem[2][0] + p04 * m.elem[3][0]);
	mat.elem[1][1] = (p00 * m.elem[0][0] + p07 * m.elem[2][0] + p08 * m.elem[3][0]) - (p01 * m.elem[0][0] + p06 * m.elem[2][0] + p09 * m.elem[3][0]);
	mat.elem[1][2] = (p03 * m.elem[0][0] + p06 * m.elem[1][0] + p11 * m.elem[3][0]) - (p02 * m.elem[0][0] + p07 * m.elem[1][0] + p10 * m.elem[3][0]);
	mat.elem[1][3] = (p04 * m.elem[0][0] + p09 * m.elem[1][0] + p10 * m.elem[2][0]) - (p05 * m.elem[0][0] + p08 * m.elem[1][0] + p11 * m.elem[2][0]);

	float q00 = m.elem[2][0] * m.elem[3][1];
	float q01 = m.elem[3][0] * m.elem[2][1];
	float q02 = m.elem[1][0] * m.elem[3][1];
	float q03 = m.elem[3][0] * m.elem[1][1];
	float q04 = m.elem[1][0] * m.elem[2][1];
	float q05 = m.elem[2][0] * m.elem[1][1];
	float q06 = m.elem[0][0] * m.elem[3][1];
	float q07 = m.elem[3][0] * m.elem[0][1];
	float q08 = m.elem[0][0] * m.elem[2][1];
	float q09 = m.elem[2][0] * m.elem[0][1];
	float q10 = m.elem[0][0] * m.elem[1][1];
	float q11 = m.elem[1][0] * m.elem[0][1];

	mat.elem[2][0] = (q00 * m.elem[1][3] + q03 * m.elem[2][3] + q04 * m.elem[3][3]) - (q01 * m.elem[1][3] + q02 * m.elem[2][3] + q05 * m.elem[3][3]);
	mat.elem[2][1] = (q01 * m.elem[0][3] + q06 * m.elem[2][3] + q09 * m.elem[3][3]) - (q00 * m.elem[0][3] + q07 * m.elem[2][3] + q08 * m.elem[3][3]);
	mat.elem[2][2] = (q02 * m.elem[0][3] + q07 * m.elem[1][3] + q10 * m.elem[3][3]) - (q03 * m.elem[0][3] + q06 * m.elem[1][3] + q11 * m.elem[3][3]);
	mat.elem[2][3] = (q05 * m.elem[0][3] + q08 * m.elem[1][3] + q11 * m.elem[2][3]) - (q04 * m.elem[0][3] + q09 * m.elem[1][3] + q10 * m.elem[2][3]);
	mat.elem[3][0] = (q02 * m.elem[2][2] + q05 * m.elem[3][2] + q01 * m.elem[1][2]) - (q04 * m.elem[3][2] + q00 * m.elem[1][2] + q03 * m.elem[2][2]);
	mat.elem[3][1] = (q08 * m.elem[3][2] + q00 * m.elem[0][2] + q07 * m.elem[2][2]) - (q06 * m.elem[2][2] + q09 * m.elem[3][2] + q01 * m.elem[0][2]);
	mat.elem[3][2] = (q06 * m.elem[1][2] + q11 * m.elem[3][2] + q03 * m.elem[0][2]) - (q10 * m.elem[3][2] + q02 * m.elem[0][2] + q07 * m.elem[1][2]);
	mat.elem[3][3] = (q10 * m.elem[2][2] + q04 * m.elem[0][2] + q09 * m.elem[1][2]) - (q08 * m.elem[1][2] + q11 * m.elem[2][2] + q05 * m.elem[0][2]);

	return mat * (1.0f / (m.elem[0][0] * mat.elem[0][0] + m.elem[1][0] * mat.elem[0][1] + m.elem[2][0] * mat.elem[0][2] + m.elem[3][0] * mat.elem[0][3]));
}


float2x2 rotate2(const float angle)
{
	float cosA = cosf(angle), sinA = sinf(angle);

	return float2x2(
		cosA, -sinA,
		sinA,  cosA);
}

float3x3 rotateX3(const float angle)
{
	float cosA = cosf(angle), sinA = sinf(angle);
	
	return float3x3(
		1, 0,     0,
		0, cosA, -sinA,
		0, sinA,  cosA);
}

float3x3 rotateY3(const float angle)
{
	float cosA = cosf(angle), sinA = sinf(angle);

	return float3x3(
		cosA, 0, -sinA,
		0,    1,  0,
		sinA, 0,  cosA);
}

float3x3 rotateZ3(const float angle)
{
	float cosA = cosf(angle), sinA = sinf(angle);

	return float3x3(
		cosA, -sinA, 0,
		sinA,  cosA, 0,
		0,     0,    1);
}

float3x3 rotateXY3(const float angleX, const float angleY)
{
	float cosX = cosf(angleX), sinX = sinf(angleX), 
		  cosY = cosf(angleY), sinY = sinf(angleY);

	return float3x3(
		 cosY,        0,    -sinY,
		-sinX * sinY, cosX, -sinX * cosY,
		 cosX * sinY, sinX,  cosX * cosY);
}

float3x3 rotateZXY3(const float angleX, const float angleY, const float angleZ)
{
	float cosX = cosf(angleX), sinX = sinf(angleX), 
		  cosY = cosf(angleY), sinY = sinf(angleY),
		  cosZ = cosf(angleZ), sinZ = sinf(angleZ);

	return float3x3(
		cosY * cosZ + sinX * sinY * sinZ,   -cosX * sinZ,    sinX * cosY * sinZ - sinY * cosZ,
		cosY * sinZ - sinX * sinY * cosZ,    cosX * cosZ,   -sinY * sinZ - sinX * cosY * cosZ,
		cosX * sinY,                         sinX,           cosX * cosY);
}

float4x4 rotateX4(const float angle)
{
	float cosA = cosf(angle), sinA = sinf(angle);
	
	return float4x4(
		1, 0,     0,    0,
		0, cosA, -sinA, 0,
		0, sinA,  cosA, 0,
		0, 0,     0,    1);
}

float4x4 rotateY4(const float angle)
{
	float cosA = cosf(angle), sinA = sinf(angle);

	return float4x4(
		cosA, 0, -sinA, 0,
		0,    1,  0,    0,
		sinA, 0,  cosA, 0,
		0,    0,  0,    1);
}

float4x4 rotateZ4(const float angle)
{
	float cosA = cosf(angle), sinA = sinf(angle);

	return float4x4(
		cosA, -sinA, 0, 0,
		sinA,  cosA, 0, 0,
		0,     0,    1, 0,
		0,     0,    0, 1);
}

float4x4 rotateXY4(const float angleX, const float angleY)
{
	float cosX = cosf(angleX), sinX = sinf(angleX), 
		  cosY = cosf(angleY), sinY = sinf(angleY);

	return float4x4(
		 cosY,        0,    -sinY,        0,
		-sinX * sinY, cosX, -sinX * cosY, 0,
		 cosX * sinY, sinX,  cosX * cosY, 0,
		 0,           0,     0,           1);
}

float4x4 rotateZXY4(const float angleX, const float angleY, const float angleZ)
{
	float cosX = cosf(angleX), sinX = sinf(angleX), 
		  cosY = cosf(angleY), sinY = sinf(angleY),
		  cosZ = cosf(angleZ), sinZ = sinf(angleZ);

	return float4x4(
		cosY * cosZ + sinX * sinY * sinZ,   -cosX * sinZ,    sinX * cosY * sinZ - sinY * cosZ,  0,
		cosY * sinZ - sinX * sinY * cosZ,    cosX * cosZ,   -sinY * sinZ - sinX * cosY * cosZ,  0,
		cosX * sinY,                         sinX,           cosX * cosY,                       0,
		0,                                   0,              0,                                 1);
}

float4x4 translate(const float x, const float y, const float z)
{
	return float4x4(1,0,0,x, 0,1,0,y, 0,0,1,z, 0,0,0,1);
}

float4x4 translate(const float3 &v)
{
	return float4x4(1,0,0,v.x, 0,1,0,v.y, 0,0,1,v.z, 0,0,0,1);
}

float2x2 scale2(const float x, const float y)
{
    return float2x2(x,0, 0,y);
}

float3x3 scale3(const float x, const float y, const float z)
{
    return float3x3(x,0,0, 0,y,0, 0,0,z);
}

float4x4 scale4(const float x, const float y, const float z)
{
    return float4x4(x,0,0,0, 0,y,0,0, 0,0,z,0, 0,0,0,1);
}

float4x4 OrthoMatrix(const float left, const float right, const float top, const float bottom, const float zNear, const float zFar)
{
	float rl = right - left;
	float tb = top - bottom;
	float fn = zFar - zNear;

	float4x4 mat(
		2.0f / rl, 0,         0,         -(right + left) / rl,
		0,         2.0f / tb, 0,         -(top + bottom) / tb,
		0,         0,        -2.0f / fn, -(zFar + zNear) / fn,
		0,         0,         0,         1);

	mat.elem[2][2] = 0.5f * (mat.elem[2][2] + mat.elem[3][2]);
	mat.elem[2][3] = 0.5f * (mat.elem[2][3] + mat.elem[3][3]);

	return mat;
}

float4x4 PerspectiveMatrixX(const float fov, const int width, const int height, const float zNear, const float zFar)
{
	float w = tanf(0.5f * fov);
	float h = (w * height) / width;

	float4x4 mat(
		1.0f / w, 0,        0, 0,
		0,        1.0f / h, 0, 0,
		0,        0,        (zFar + zNear) / (zFar - zNear), -(2 * zFar * zNear) / (zFar - zNear),
		0,        0,        1, 0);

	mat.elem[2][2] = 0.5f * (mat.elem[2][2] + mat.elem[3][2]);
	mat.elem[2][3] = 0.5f * (mat.elem[2][3] + mat.elem[3][3]);

	return mat;
}

float4x4 PerspectiveMatrixY(const float fov, const int width, const int height, const float zNear, const float zFar)
{
	float h = tanf(0.5f * fov);
	float w = (h * width) / height;

	float4x4 mat(
		1.0f / w, 0,        0, 0,
		0,        1.0f / h, 0, 0,
		0,        0,        (zFar + zNear) / (zFar - zNear), -(2 * zFar * zNear) / (zFar - zNear),
		0,        0,        1, 0);

	mat.elem[2][2] = 0.5f * (mat.elem[2][2] + mat.elem[3][2]);
	mat.elem[2][3] = 0.5f * (mat.elem[2][3] + mat.elem[3][3]);

	return mat;
}

float4x4 ShadowMatrix(const float3 &planeNormal, const float planeOffset, const float3 &lightPos)
{
	float dist = dot(lightPos, planeNormal) + planeOffset;

	return float4x4(
		dist - lightPos.x * planeNormal.x,
		     - lightPos.x * planeNormal.y,
		     - lightPos.x * planeNormal.z,
		     - lightPos.x * planeOffset,

		     - lightPos.y * planeNormal.x,
		dist - lightPos.y * planeNormal.y,
		     - lightPos.y * planeNormal.z,
		     - lightPos.y * planeOffset,

		     - lightPos.z * planeNormal.x,
		     - lightPos.z * planeNormal.y,
		dist - lightPos.z * planeNormal.z,
		     - lightPos.z * planeOffset,

		     - planeNormal.x,
		     - planeNormal.y,
		     - planeNormal.z,
		dist - planeOffset);
}

float4x4 MirrorMatrix(const float3 &planeNormal, const float planeOffset)
{
	return float4x4(
		1 - 2 * planeNormal.x * planeNormal.x,
		  - 2 * planeNormal.x * planeNormal.y,
		  - 2 * planeNormal.x * planeNormal.z,
		  - 2 * planeNormal.x * planeOffset,
		
		  - 2 * planeNormal.y * planeNormal.x,
		1 - 2 * planeNormal.y * planeNormal.y,
		  - 2 * planeNormal.y * planeNormal.z,
		  - 2 * planeNormal.y * planeOffset,

		  - 2 * planeNormal.z * planeNormal.x,
		  - 2 * planeNormal.z * planeNormal.y,
		1 - 2 * planeNormal.z * planeNormal.z,
		  - 2 * planeNormal.z * planeOffset,

		0, 0, 0, 1);
}

float4x4 CubemapModelviewMatrix(const unsigned int face)
{
	switch (face)
	{
		case 0:
			return float4x4(
				0, 0, -1, 0,
				0, 1,  0, 0,
				1, 0,  0, 0,
				0, 0,  0, 1);
		case 1:
			return float4x4(
				 0, 0, 1, 0,
				 0, 1, 0, 0,
				-1, 0, 0, 0,
				 0, 0, 0, 1);
		case 2:
			return float4x4(
				1, 0,  0, 0,
				0, 0, -1, 0,
				0, 1,  0, 0,
				0, 0,  0, 1);
		case 3:
			return float4x4(
				1,  0, 0, 0,
				0,  0, 1, 0,
				0, -1, 0, 0,
				0,  0, 0, 1);
		case 4:
			return float4x4(
				1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1);
		default:
			return float4x4(
				-1, 0,  0, 0,
				 0, 1,  0, 0,
				 0, 0, -1, 0,
				 0, 0,  0, 1);
	}
}

float4x4 CubemapProjectionMatrix(const float zNear, const float zFar)
{
	float4x4 mat(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, (zFar + zNear) / (zFar - zNear), -(2 * zFar * zNear) / (zFar - zNear),
		0, 0, 1, 0);

	mat.elem[2][2] = 0.5f * (mat.elem[2][2] + mat.elem[3][2]);
	mat.elem[2][3] = 0.5f * (mat.elem[2][3] + mat.elem[3][3]);

	return mat;
}

float4x4 LookAt(const float3 &from, const float3 &at)
{
	float3 z = normalize(at - from);
	float3 x = normalize(float3(z.z, 0, -z.x));
	float3 y = normalize(cross(z, x));

	float4x4 mat(
		x.x, x.y, x.z, 0,
		y.x, y.y, y.z, 0,
		z.x, z.y, z.z, 0,
		0,   0,   0,   1);
	mat.translate(-from);

	return mat;
}

float4x4 OpenGLToD3DProjectionMatrix(const float4x4 &m)
{
	float4x4 mat = m;

	mat.elem[2][0] = 0.5f * (mat.elem[2][0] + mat.elem[3][0]);
	mat.elem[2][1] = 0.5f * (mat.elem[2][1] + mat.elem[3][1]);
	mat.elem[2][2] = 0.5f * (mat.elem[2][2] + mat.elem[3][2]);
	mat.elem[2][3] = 0.5f * (mat.elem[2][3] + mat.elem[3][3]);

	return mat;
}

float4x4 PegToFarPlane(const float4x4 &m)
{
	float4x4 mat;

	memcpy(mat.elem[0], m.elem[0], 8 * sizeof(float));
	memcpy(mat.elem[2], m.elem[3], 4 * sizeof(float));
	memcpy(mat.elem[3], m.elem[3], 4 * sizeof(float));

	return mat;
}

float2x2 Identity2()
{
	return float2x2(1,0, 0,1);
}

float3x3 Identity3()
{
	return float3x3(1,0,0, 0,1,0, 0,0,1);
}

float4x4 Identity4()
{
	return float4x4(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
}
