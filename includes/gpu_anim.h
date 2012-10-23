/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and
 * proprietary rights in and to this software and related documentation.
 * Any use, reproduction, disclosure, or distribution of this software
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA)
 * associated with this source code for terms and conditions that govern
 * your use of this NVIDIA software.
 *
 */


#ifndef __GPU_ANIM_H__
#define __GPU_ANIM_H__

#include "gl_helper.h"

#include "cuda_by_example.h"
#include "cuda.h"
#include "cuda_gl_interop.h"
#include <iostream>
#include "cutil_math.h"
#include "../src/main.cuh"
#include "../src/functions.cuh"


PFNGLBINDBUFFERARBPROC    glBindBuffer     = NULL;
PFNGLDELETEBUFFERSARBPROC glDeleteBuffers  = NULL;
PFNGLGENBUFFERSARBPROC    glGenBuffers     = NULL;
PFNGLBUFFERDATAARBPROC    glBufferData     = NULL;


struct GPUAnimBitmap {
    GLuint  bufferObj;
    cudaGraphicsResource *resource;
    int     width, height;
    void    *dataBlock;
    void (*fAnim)(uchar4*,void*,int);
    void (*animExit)(void*);
    void (*clickDrag)(void*,int,int,int,int);
    int     dragStartX, dragStartY;

    GPUAnimBitmap( int w, int h, void *d = NULL ) {
        width = w;
        height = h;
        dataBlock = d;
        clickDrag = NULL;

        // first, find a CUDA device and set it to graphic interop
        cudaDeviceProp  prop;
        int dev;
        memset( &prop, 0, sizeof( cudaDeviceProp ) );
        prop.major = 1;
        prop.minor = 0;
        HANDLE_ERROR( cudaChooseDevice( &dev, &prop ) );
        cudaGLSetGLDevice( dev );

        // a bug in the Windows GLUT implementation prevents us from
        // passing zero arguments to glutInit()
        int c=1;
        char* dummy = "";
        glutInit( &c, &dummy );
        glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA );
        glutInitWindowSize( width, height );
        glutCreateWindow( "bitmap" );

        glBindBuffer    = (PFNGLBINDBUFFERARBPROC)GET_PROC_ADDRESS("glBindBuffer");
        glDeleteBuffers = (PFNGLDELETEBUFFERSARBPROC)GET_PROC_ADDRESS("glDeleteBuffers");
        glGenBuffers    = (PFNGLGENBUFFERSARBPROC)GET_PROC_ADDRESS("glGenBuffers");
        glBufferData    = (PFNGLBUFFERDATAARBPROC)GET_PROC_ADDRESS("glBufferData");

        glGenBuffers( 1, &bufferObj );
        glBindBuffer( GL_PIXEL_UNPACK_BUFFER_ARB, bufferObj );
        glBufferData( GL_PIXEL_UNPACK_BUFFER_ARB, width * height * 4,
                      NULL, GL_DYNAMIC_DRAW_ARB );

        HANDLE_ERROR( cudaGraphicsGLRegisterBuffer( &resource, bufferObj, cudaGraphicsMapFlagsNone ) );
    }

    ~GPUAnimBitmap() {
        free_resources();
    }

    void free_resources( void ) {
        HANDLE_ERROR( cudaGraphicsUnregisterResource( resource ) );

        glBindBuffer( GL_PIXEL_UNPACK_BUFFER_ARB, 0 );
        glDeleteBuffers( 1, &bufferObj );
    }


    long image_size( void ) const { return width * height * 4; }

    void click_drag( void (*f)(void*,int,int,int,int)) {
        clickDrag = f;
    }

    void anim_and_exit( void (*f)(uchar4*,void*,int), void(*e)(void*) ) {
        GPUAnimBitmap**   bitmap = get_bitmap_ptr();
        *bitmap = this;
        fAnim = f;
        animExit = e;

        glutKeyboardFunc( Key );
        glutDisplayFunc( Draw );
        if (clickDrag != NULL)
            glutMouseFunc( mouse_func );
        glutIdleFunc( idle_func );
        glutMainLoop();
        
    }

    // static method used for glut callbacks
    static GPUAnimBitmap** get_bitmap_ptr( void ) {
        static GPUAnimBitmap*   gBitmap;
        return &gBitmap;
    }

    // static method used for glut callbacks
    static void mouse_func( int button, int state,
                            int mx, int my ) {
        if (button == GLUT_LEFT_BUTTON) {
            GPUAnimBitmap*   bitmap = *(get_bitmap_ptr());
            if (state == GLUT_DOWN) {
                bitmap->dragStartX = mx;
                bitmap->dragStartY = my;
            } else if (state == GLUT_UP) {
                bitmap->clickDrag( bitmap->dataBlock,
                                   bitmap->dragStartX,
                                   bitmap->dragStartY,
                                   mx, my );
            }
        }
    }

    // static method used for glut callbacks
    static void idle_func( void ) {
        static int ticks = 1;
        GPUAnimBitmap*  bitmap = *(get_bitmap_ptr());
        uchar4*         devPtr;
        size_t  size;
		DataBlock*		simBlock = (DataBlock*) bitmap->dataBlock;

        HANDLE_ERROR( cudaGraphicsMapResources( 1, &(bitmap->resource), NULL ) );
        HANDLE_ERROR( cudaGraphicsResourceGetMappedPointer( (void**)&devPtr, &size, bitmap->resource) );

        bitmap->fAnim( devPtr, bitmap->dataBlock, ticks++ );

        HANDLE_ERROR( cudaGraphicsUnmapResources( 1, &(bitmap->resource), NULL ) );

        glutPostRedisplay();
        
        // Testando se o tempo de simulação já foi atingido (-1 = tempo infinito)
		if (simBlock->timeCtrl.simDuration != -1 && simBlock->timeCtrl.tempo >= simBlock->timeCtrl.simDuration) Key(27,0,0);
    }

    // static method used for glut callbacks
    static void Key(unsigned char key, int x, int y)
	{
		GPUAnimBitmap*		bitmap = *(get_bitmap_ptr());
		DataBlock*			simBlock = (DataBlock*) bitmap->dataBlock;
		SystemProperties*	sisProps = &simBlock->sisProps;
		RenderParameters*	renderPar = &simBlock->renderPar;

		int*	bgColor = &renderPar->bgColor;
		float2*	gravity = &sisProps->gravity;
		bool*	viewRotations = &renderPar->viewRotations;
		bool*	colorByPressure = &renderPar->colorByPressure;
		float*	pressure = &simBlock->partValues.pressure[0];
		int		numParts = sisProps->numParticles;

        switch (key)
		{
			case 'b':
				if (*bgColor == 0) *bgColor = 255;
				else *bgColor = 0;
				break;

			case 'g':
				if (gravity->y != 0) set_gravity (sisProps, make_float2(0.0f));
				else set_gravity (sisProps, make_float2(0.0, -9.81f));
				break;

			case 'r':
				if (*viewRotations) set_viewRotations (renderPar, false);
				else set_viewRotations (renderPar, true);
				break;

			case 'p':
				if (*colorByPressure) set_colorByPressure (renderPar, false);
				else set_colorByPressure (renderPar, true);
				break;

			case 'u':
				updatePressureScale (renderPar, pressure, numParts);
				break;

            case 27:
                if (bitmap->animExit)
                    bitmap->animExit( bitmap->dataBlock );
                bitmap->free_resources();
                exit(0);
        }
    }

    // static method used for glut callbacks
    static void Draw( void ) {
        GPUAnimBitmap*   bitmap = *(get_bitmap_ptr());
        glClearColor( 0.0, 0.0, 0.0, 1.0 );
        glClear( GL_COLOR_BUFFER_BIT );
        glDrawPixels( bitmap->width, bitmap->height, GL_RGBA,
                      GL_UNSIGNED_BYTE, 0 );
        glutSwapBuffers();
    }
};


#endif  // __GPU_ANIM_H__

