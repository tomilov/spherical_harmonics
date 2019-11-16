//-----------------------------------------------------------------------------
// File: Framework\Platform.h
// Copyright (c) 2007 Advanced Micro Devices Inc. All rights reserved.
//-----------------------------------------------------------------------------





#ifndef _PLATFORM_H_
#define _PLATFORM_H_

#include <stddef.h>
#include <tchar.h>

#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

/** \file
    Platform definitions.  This file defines various basic types and a handful use useful macros.
*/

// To get rid of MSVC 2005 warnings
#if _MSC_VER >= 1400
#define _CRT_SECURE_NO_DEPRECATE
#pragma warning (disable: 4996)
#endif

/** \defgroup paths Media Path Definitions
    Provides a simple predefined path for media files.
    @{
*/
#define MEDIA_PATH _T("../../Media/")
#define TEX_PATH    MEDIA_PATH _T("Textures/")
#define MODEL_PATH  MEDIA_PATH _T("Models/")
#define FONT_PATH   MEDIA_PATH _T("Fonts/")
#define SHADER_PATH MEDIA_PATH _T("Shaders/")
#define MUSIC_PATH  MEDIA_PATH _T("Music/")
//@}

/** \defgroup types Data Types
    Core data type definition.
    @{
*/
typedef int int16;               /**< signed 16-bit integer */
typedef int int32;               /**< signed 32-bit integer */
typedef unsigned short uint16;   /**< unsigned 16-bit integer */
typedef unsigned int uint32;     /**< unsigned 32-bit integer */
typedef unsigned int uint;       /**< unsigned integer (same as uint32 on 32-bit platforms) */
typedef unsigned short ushort;   /**< unsigned short */
typedef unsigned char ubyte;     /**< unsigned byte */
typedef __int64 int64;           /**< signed 64-bit integer */
typedef unsigned __int64 uint64; /**< unsigned 64-bit integer */
typedef intptr_t intptr;         /**< pointer to a signed integer */
typedef uintptr_t uintptr;       /**< pointer to an unsigned integer */
//@}

/** 16-bit float represented as an unsigned short */
struct half
{
	unsigned short sh;

	half(){};
	half(const float x);
	operator float () const;
};



#define PI 3.1415926535897932f


/** \defgroup platformUtilMacro Useful Macros
    Just some handy utility macros.
    @{
*/
#define elementsOf(x) (sizeof(x) / sizeof(x[0]))                            /**< returns the number of elements in an array of items */
#define offsetOf(strct, x) intptr(&((strct *) NULL)->x)                     /**< returns a pointer to an element in a block of memory */
#define STR(x) #x
#define DEF_STR(x) STR(x)
#define DEF_MACRO(x) { #x, DEF_STR(x) }
//@}
#ifndef NULL
#define NULL 0
#endif

#if defined(_DEBUG) && !defined(DEBUG)
#define DEBUG
#endif

#pragma warning(push)
#pragma warning(disable: 4035)
inline uint64 GetCycle(){
    return __rdtsc();
}
#pragma warning(pop)

#endif // _PLATFORM_H_
