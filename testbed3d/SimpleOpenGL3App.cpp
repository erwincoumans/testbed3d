#include "SimpleOpenGL3App.h"


#ifdef __APPLE__
#include "OpenGLWindow/MacOpenGLWindow.h"
#else

#include "GL/glew.h"
#ifdef _WIN32
#include "OpenGLWindow/Win32OpenGLWindow.h"
#else
//let's cross the fingers it is Linux/X11
#include "OpenGLWindow/X11OpenGLWindow.h"
#endif //_WIN32
#endif//__APPLE__

#include "OpenGLWindow/GLPrimitiveRenderer.h"
#include "OpenGLWindow/GLInstancingRenderer.h"

#include "Bullet3Common/b3Vector3.h"

#include "../btgui/OpenGLTrueTypeFont/fontstash.h"
#include "../btgui/OpenGLWindow/TwFonts.h"
#include "../btgui/Bullet3Common/b3Transform.h"

struct SimpleInternalData
{
	GLuint m_fontTextureId; 
};

static SimpleOpenGL3App* gApp=0;

void SimpleResizeCallback( float width, float height)
{
	gApp->m_instancingRenderer->resize(width,height);
	gApp->m_primRenderer->setScreenSize(width,height);
	
}

static GLuint BindFont(const CTexFont *_Font)
{
    GLuint TexID = 0;
    glGenTextures(1, &TexID);
    glBindTexture(GL_TEXTURE_2D, TexID);
    glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, _Font->m_TexWidth, _Font->m_TexHeight, 0, GL_RED, GL_UNSIGNED_BYTE, _Font->m_TexBytes);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);

    return TexID;
}


SimpleOpenGL3App::SimpleOpenGL3App(	const char* title, int width,int height)
{
	gApp = this;
	m_data = new SimpleInternalData;
	m_window = new b3gDefaultOpenGLWindow();
	b3gWindowConstructionInfo ci;
	ci.m_title = title;
	ci.m_width = width;
	ci.m_height = height;
	m_window->createWindow(ci);
	
	m_window->setWindowTitle(title);
	glClearColor(1,1,1,1);
	m_window->startRendering();
#ifndef __APPLE__
	glewInit();
#endif

	m_primRenderer = new GLPrimitiveRenderer(width,height);
	
	m_instancingRenderer = new GLInstancingRenderer(128*1024,4*1024*1024);
	m_instancingRenderer->init();
	m_instancingRenderer->resize(width,height);
	m_instancingRenderer->InitShaders();

	m_window->setMouseMoveCallback(b3DefaultMouseMoveCallback);
	m_window->setMouseButtonCallback(b3DefaultMouseButtonCallback);
	m_window->setKeyboardCallback(b3DefaultKeyboardCallback);
	m_window->setWheelCallback(b3DefaultWheelCallback);
	m_window->setResizeCallback(SimpleResizeCallback);

	TwGenerateDefaultFonts();
	m_data->m_fontTextureId = BindFont(g_DefaultNormalFont);
}


void	SimpleOpenGL3App::drawLine(const float* from,const float* to)
{
	
	int numV = m_vertices.size();
	m_vertices.push_back(b3MakeVector3(from[0],from[1],from[2]));
	m_vertices.push_back(b3MakeVector3(to[0],to[1],to[2]));
	m_indices.push_back(numV);
	m_indices.push_back(numV+1);
}




void SimpleOpenGL3App::drawArc(const b3Vector3& center, const b3Vector3& normal, const b3Vector3& axis, b3Scalar radiusA, b3Scalar radiusB, b3Scalar minAngle, b3Scalar maxAngle, 
			 bool drawSect, b3Scalar stepDegrees )
{
	const b3Vector3& vx = axis;
	b3Vector3 vy = normal.cross(axis);
	b3Scalar step = stepDegrees * B3_RADS_PER_DEG;
	int nSteps = (int)((maxAngle - minAngle) / step);
	if(!nSteps) nSteps = 1;
	b3Vector3 prev = center + radiusA * vx * b3Cos(minAngle) + radiusB * vy * b3Sin(minAngle);
	if(drawSect)
	{
		drawLine(center, prev);
	}
	for(int i = 1; i <= nSteps; i++)
	{
		b3Scalar angle = minAngle + (maxAngle - minAngle) * b3Scalar(i) / b3Scalar(nSteps);
		b3Vector3 next = center + radiusA * vx * b3Cos(angle) + radiusB * vy * b3Sin(angle);
		drawLine(prev, next);
		prev = next;
	}
	if(drawSect)
	{
		drawLine(center, prev);
	}
}
	

void SimpleOpenGL3App::drawCylinder(float radius, float halfHeight, int upAxis, const float* pos, const float* ornQuat)
{
	b3Transform transform;
	transform.setIdentity();
	transform.setOrigin(b3MakeVector3(pos[0],pos[1],pos[2]));
	transform.setRotation(b3Quaternion(ornQuat[0],ornQuat[1],ornQuat[2],ornQuat[3]));

	b3Vector3 start = transform.getOrigin();
	b3Vector3	offsetHeight=b3MakeVector3(0,0,0);
	offsetHeight[upAxis] = halfHeight;
	int stepDegrees=30;
	b3Vector3 capStart=b3MakeVector3(0.f,0.f,0.f);
	capStart[upAxis] = -halfHeight;
	b3Vector3 capEnd=b3MakeVector3(0.f,0.f,0.f);
	capEnd[upAxis] = halfHeight;

	for (int i=0;i<360;i+=stepDegrees)
	{
		capEnd[(upAxis+1)%3] = capStart[(upAxis+1)%3] = b3Sin(b3Scalar(i)*B3_RADS_PER_DEG)*radius;
		capEnd[(upAxis+2)%3] = capStart[(upAxis+2)%3]  = b3Cos(b3Scalar(i)*B3_RADS_PER_DEG)*radius;
		drawLine(start+transform.getBasis() * capStart,start+transform.getBasis() * capEnd);
	}
	// Drawing top and bottom caps of the cylinder
	b3Vector3 yaxis=b3MakeVector3(0,0,0);
	yaxis[upAxis] = b3Scalar(1.0);
	b3Vector3 xaxis=b3MakeVector3(0,0,0);
	xaxis[(upAxis+1)%3] = b3Scalar(1.0);
	drawArc(start-transform.getBasis()*(offsetHeight),transform.getBasis()*yaxis,transform.getBasis()*xaxis,radius,radius,0,B3_2_PI,false,b3Scalar(10.0));
	drawArc(start+transform.getBasis()*(offsetHeight),transform.getBasis()*yaxis,transform.getBasis()*xaxis,radius,radius,0,B3_2_PI,false,b3Scalar(10.0));
}

void SimpleOpenGL3App::drawCone(float radius, float height, int upAxis, const float* pos, const float* ornQuat)
{
	b3Transform transform;
	transform.setIdentity();
	transform.setOrigin(b3MakeVector3(pos[0],pos[1],pos[2]));
	transform.setRotation(b3Quaternion(ornQuat[0],ornQuat[1],ornQuat[2],ornQuat[3]));

	int stepDegrees = 30;
	b3Vector3 start = transform.getOrigin();

	b3Vector3	offsetHeight=b3MakeVector3(0,0,0);
	b3Scalar halfHeight = height * b3Scalar(0.5);
	offsetHeight[upAxis] = halfHeight;
	b3Vector3	offsetRadius=b3MakeVector3(0,0,0);
	offsetRadius[(upAxis+1)%3] = radius;
	b3Vector3	offset2Radius=b3MakeVector3(0,0,0);
	offset2Radius[(upAxis+2)%3] = radius;


	b3Vector3 capEnd=b3MakeVector3(0.f,0.f,0.f);
	capEnd[upAxis] = -halfHeight;

	for (int i=0;i<360;i+=stepDegrees)
	{
		capEnd[(upAxis+1)%3] = b3Sin(b3Scalar(i)*B3_RADS_PER_DEG)*radius;
		capEnd[(upAxis+2)%3] = b3Cos(b3Scalar(i)*B3_RADS_PER_DEG)*radius;
		drawLine(start+transform.getBasis() * (offsetHeight),start+transform.getBasis() * capEnd);
	}

	drawLine(start+transform.getBasis() * (offsetHeight),start+transform.getBasis() * (-offsetHeight+offsetRadius));
	drawLine(start+transform.getBasis() * (offsetHeight),start+transform.getBasis() * (-offsetHeight-offsetRadius));
	drawLine(start+transform.getBasis() * (offsetHeight),start+transform.getBasis() * (-offsetHeight+offset2Radius));
	drawLine(start+transform.getBasis() * (offsetHeight),start+transform.getBasis() * (-offsetHeight-offset2Radius));

	// Drawing the base of the cone
	b3Vector3 yaxis=b3MakeVector3(0,0,0);
	yaxis[upAxis] = b3Scalar(1.0);
	b3Vector3 xaxis=b3MakeVector3(0,0,0);
	xaxis[(upAxis+1)%3] = b3Scalar(1.0);
	drawArc(start-transform.getBasis()*(offsetHeight),transform.getBasis()*yaxis,transform.getBasis()*xaxis,radius,radius,0,B3_2_PI,false,10.0);
}

void SimpleOpenGL3App::drawSpherePatch(const b3Vector3& center, const b3Vector3& up, const b3Vector3& axis, b3Scalar radius, 
		b3Scalar minTh, b3Scalar maxTh, b3Scalar minPs, b3Scalar maxPs, b3Scalar stepDegrees ,bool drawCenter)
	{
		b3Vector3 vA[74];
		b3Vector3 vB[74];
		b3Vector3 *pvA = vA, *pvB = vB, *pT;
		b3Vector3 npole = center + up * radius;
		b3Vector3 spole = center - up * radius;
		b3Vector3 arcStart;
		b3Scalar step = stepDegrees * B3_RADS_PER_DEG;
		const b3Vector3& kv = up;
		const b3Vector3& iv = axis;
		b3Vector3 jv = kv.cross(iv);
		bool drawN = false;
		bool drawS = false;
		if(minTh <= -B3_HALF_PI)
		{
			minTh = -B3_HALF_PI + step;
			drawN = true;
		}
		if(maxTh >= B3_HALF_PI)
		{
			maxTh = B3_HALF_PI - step;
			drawS = true;
		}
		if(minTh > maxTh)
		{
			minTh = -B3_HALF_PI + step;
			maxTh =  B3_HALF_PI - step;
			drawN = drawS = true;
		}
		int n_hor = (int)((maxTh - minTh) / step) + 1;
		if(n_hor < 2) n_hor = 2;
		b3Scalar step_h = (maxTh - minTh) / b3Scalar(n_hor - 1);
		bool isClosed = false;
		if(minPs > maxPs)
		{
			minPs = -B3_PI + step;
			maxPs =  B3_PI;
			isClosed = true;
		}
		else if((maxPs - minPs) >= B3_PI * b3Scalar(2.f))
		{
			isClosed = true;
		}
		else
		{
			isClosed = false;
		}
		int n_vert = (int)((maxPs - minPs) / step) + 1;
		if(n_vert < 2) n_vert = 2;
		b3Scalar step_v = (maxPs - minPs) / b3Scalar(n_vert - 1);
		for(int i = 0; i < n_hor; i++)
		{
			b3Scalar th = minTh + b3Scalar(i) * step_h;
			b3Scalar sth = radius * b3Sin(th);
			b3Scalar cth = radius * b3Cos(th);
			for(int j = 0; j < n_vert; j++)
			{
				b3Scalar psi = minPs + b3Scalar(j) * step_v;
				b3Scalar sps = b3Sin(psi);
				b3Scalar cps = b3Cos(psi);
				pvB[j] = center + cth * cps * iv + cth * sps * jv + sth * kv;
				if(i)
				{
					drawLine(pvA[j], pvB[j]);
				}
				else if(drawS)
				{
					drawLine(spole, pvB[j]);
				}
				if(j)
				{
					drawLine(pvB[j-1], pvB[j]);
				}
				else
				{
					arcStart = pvB[j];
				}
				if((i == (n_hor - 1)) && drawN)
				{
					drawLine(npole, pvB[j]);
				}
				
				if (drawCenter)
				{
					if(isClosed)
					{
						if(j == (n_vert-1))
						{
							drawLine(arcStart, pvB[j]);
						}
					}
					else
					{
						if(((!i) || (i == (n_hor-1))) && ((!j) || (j == (n_vert-1))))
						{
							drawLine(center, pvB[j]);
						}
					}
				}
			}
			pT = pvA; pvA = pvB; pvB = pT;
		}
	}

void SimpleOpenGL3App::drawCapsule(float radius, float halfHeight, int upAxis,  const float* pos, const float* ornQuat)
{
	b3Transform transform;
	transform.setIdentity();
	transform.setOrigin(b3MakeVector3(pos[0],pos[1],pos[2]));
	transform.setRotation(b3Quaternion(ornQuat[0],ornQuat[1],ornQuat[2],ornQuat[3]));

	int stepDegrees = 30;

	b3Vector3 capStart = b3MakeVector3(0.f,0.f,0.f);
	capStart[upAxis] = -halfHeight;

	b3Vector3 capEnd=b3MakeVector3(0.f,0.f,0.f);
	capEnd[upAxis] = halfHeight;

	// Draw the ends
	{

		b3Transform childTransform = transform;
		childTransform.getOrigin() = transform * capStart;
		{
			b3Vector3 center = childTransform.getOrigin();
			b3Vector3 up = childTransform.getBasis().getColumn((upAxis+1)%3);
			b3Vector3 axis = -childTransform.getBasis().getColumn(upAxis);
			b3Scalar minTh = -B3_HALF_PI;
			b3Scalar maxTh = B3_HALF_PI;
			b3Scalar minPs = -B3_HALF_PI;
			b3Scalar maxPs = B3_HALF_PI;
				
			drawSpherePatch(center, up, axis, radius,minTh, maxTh, minPs, maxPs, b3Scalar(stepDegrees) ,false);
		}



	}

	{
		b3Transform childTransform = transform;
		childTransform.getOrigin() = transform * capEnd;
		{
			b3Vector3 center = childTransform.getOrigin();
			b3Vector3 up = childTransform.getBasis().getColumn((upAxis+1)%3);
			b3Vector3 axis = childTransform.getBasis().getColumn(upAxis);
			b3Scalar minTh = -B3_HALF_PI;
			b3Scalar maxTh = B3_HALF_PI;
			b3Scalar minPs = -B3_HALF_PI;
			b3Scalar maxPs = B3_HALF_PI;
			drawSpherePatch(center, up, axis, radius,minTh, maxTh, minPs, maxPs, b3Scalar(stepDegrees) ,false);
		}
	}

	// Draw some additional lines
	b3Vector3 start = transform.getOrigin();

	for (int i=0;i<360;i+=stepDegrees)
	{
		capEnd[(upAxis+1)%3] = capStart[(upAxis+1)%3] = b3Sin(b3Scalar(i)*B3_RADS_PER_DEG)*radius;
		capEnd[(upAxis+2)%3] = capStart[(upAxis+2)%3]  = b3Cos(b3Scalar(i)*B3_RADS_PER_DEG)*radius;
		drawLine(start+transform.getBasis() * capStart,start+transform.getBasis() * capEnd);
	}
		
}

 void	SimpleOpenGL3App::drawSphere(float radius,  const float* pos, const float* ornQuat)
{
	b3Transform transform;
	transform.setIdentity();
	transform.setOrigin(b3MakeVector3(pos[0],pos[1],pos[2]));
	transform.setRotation(b3Quaternion(ornQuat[0],ornQuat[1],ornQuat[2],ornQuat[3]));

		
	b3Vector3 center = transform.getOrigin();
	b3Vector3 up = transform.getBasis().getColumn(1);
	b3Vector3 axis = transform.getBasis().getColumn(0);
	b3Scalar minTh = -B3_HALF_PI;
	b3Scalar maxTh = B3_HALF_PI;
	b3Scalar minPs = -B3_HALF_PI;
	b3Scalar maxPs = B3_HALF_PI;
	b3Scalar stepDegrees = 30.f;
	drawSpherePatch(center, up, axis, radius,minTh, maxTh, minPs, maxPs,  stepDegrees ,false);
	drawSpherePatch(center, up, -axis, radius,minTh, maxTh, minPs, maxPs, stepDegrees,false );
}


void SimpleOpenGL3App::drawBox(const float* halfDimensions, const float* position, const float* orientation)
{
	float pointColor[4]={1,0.4,0.4,1};
	//gRenderer->drawLine(b3MakeVector3(0,0,0),b3MakeVector3(1,0,0),b3MakeVector3(1,0,0),3);

	b3Vector3 org=b3MakeVector3(position[0],position[1],position[2]);
	b3Quaternion orn(orientation[0],orientation[1],orientation[2],orientation[3]);
	b3Vector3 dimx=b3MakeVector3(halfDimensions[0],0,0);
	b3Vector3 dimy=b3MakeVector3(0,halfDimensions[1],0);
	b3Vector3 dimz=b3MakeVector3(0,0,halfDimensions[2]);

	b3Transform tr(orn,org);
	b3Vector3 dx = tr.getBasis()*dimx;
	b3Vector3 dy = tr.getBasis()*dimy;
	b3Vector3 dz = tr.getBasis()*dimz;

	int startIndex=m_vertices.size();

	m_vertices.push_back(org - dx - dy - dz);
	m_vertices.push_back(org + dx - dy - dz);
	m_vertices.push_back(org + dx - dy + dz);
	m_vertices.push_back(org - dx - dy + dz);
	m_vertices.push_back(org - dx + dy - dz);
	m_vertices.push_back(org + dx + dy - dz);
	m_vertices.push_back(org + dx + dy + dz);
	m_vertices.push_back(org - dx + dy + dz);
	
	float lineColor[4]={0.,0.,.1,1};

	float lineSize = 2;
	
	
	m_indices.push_back(startIndex+0); m_indices.push_back(startIndex+1);
	m_indices.push_back(startIndex+1); m_indices.push_back(startIndex+2);
	m_indices.push_back(startIndex+2); m_indices.push_back(startIndex+3);
	m_indices.push_back(startIndex+3); m_indices.push_back(startIndex+0);

	m_indices.push_back(startIndex+4); m_indices.push_back(startIndex+5);
	m_indices.push_back(startIndex+5); m_indices.push_back(startIndex+6);
	m_indices.push_back(startIndex+6); m_indices.push_back(startIndex+7);
	m_indices.push_back(startIndex+7); m_indices.push_back(startIndex+4);

	m_indices.push_back(startIndex+0); m_indices.push_back(startIndex+4);
	m_indices.push_back(startIndex+1); m_indices.push_back(startIndex+5);
	m_indices.push_back(startIndex+2); m_indices.push_back(startIndex+6);
	m_indices.push_back(startIndex+3); m_indices.push_back(startIndex+7);

	

}


void SimpleOpenGL3App::drawText( const char* txt, int posX, int posY)
{
		
    
    
        
    //
    //printf("str = %s\n",unicodeText);
    int xpos=0;
    int ypos=0;
    float dx;
        
    int measureOnly=0;
        
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	/*
	if (m_useTrueTypeFont)
	{
			
		float yoffset = 0.f;
		if (m_retinaScale==2.0f)
		{
			yoffset = -12;
		}
		Translate(r);
		sth_draw_text(m_font,
                    1,m_fontScaling,
                    r.x,r.y+yoffset,
                    unicodeText,&dx, m_screenWidth,m_screenHeight,measureOnly,m_retinaScale);
			 
	} else
	*/
	{
		//float width = 0.f;
		int pos=0;
		float color[]={0.2f,0.2,0.2f,1.f};
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D,m_data->m_fontTextureId);
		
		//float width = r.x;
		float extraSpacing = 0.;

		int startX = posX;
		int startY = posY;
		

		while (txt[pos])
		{
			int c = txt[pos];
			//r.h = g_DefaultNormalFont->m_CharHeight;
			//r.w = g_DefaultNormalFont->m_CharWidth[c]+extraSpacing;
			int endX = startX+g_DefaultNormalFont->m_CharWidth[c];
			int endY = startY+g_DefaultNormalFont->m_CharHeight;
			//Gwen::Rect rect = r;
			//Translate( rect );

			
			float currentColor[]={0.2f,0.2,0.2f,1.f};

			m_primRenderer->drawTexturedRect(startX, startY, endX, endY, currentColor,g_DefaultNormalFont->m_CharU0[c],g_DefaultNormalFont->m_CharV0[c],g_DefaultNormalFont->m_CharU1[c],g_DefaultNormalFont->m_CharV1[c]);

			//DrawTexturedRect(0,r,g_DefaultNormalFont->m_CharU0[c],g_DefaultNormalFont->m_CharV0[c],g_DefaultNormalFont->m_CharU1[c],g_DefaultNormalFont->m_CharV1[c]);
		//	DrawFilledRect(r);

			startX = endX;
			//startY = endY;
			
			pos++;
				
		}
		glBindTexture(GL_TEXTURE_2D,0);
	}

	glDisable(GL_BLEND);
}

void SimpleOpenGL3App::drawGrid(int gridSize)
{
	
	b3Vector3 gridColor = b3MakeVector3(0.5,0.5,0.5);
	for(int i=-gridSize;i<=gridSize;i++)
	{
		
		GLint err = glGetError();
		b3Assert(err==GL_NO_ERROR);
		
		m_instancingRenderer->drawLine(b3MakeVector3(float(i),0,float(-gridSize)),b3MakeVector3(float(i),0,float(gridSize)),gridColor);
		
		err = glGetError();
		b3Assert(err==GL_NO_ERROR);
		
		m_instancingRenderer->drawLine(b3MakeVector3(float(-gridSize),0,float(i)),b3MakeVector3(float(gridSize),0,float(i)),gridColor);
	}
				
	m_instancingRenderer->drawLine(b3MakeVector3(0,0,0),b3MakeVector3(1,0,0),b3MakeVector3(1,0,0),3);
	m_instancingRenderer->drawLine(b3MakeVector3(0,0,0),b3MakeVector3(0,1,0),b3MakeVector3(0,1,0),3);
	m_instancingRenderer->drawLine(b3MakeVector3(0,0,0),b3MakeVector3(0,0,1),b3MakeVector3(0,0,1),3);

	m_instancingRenderer->drawPoint(b3MakeVector3(1,0,0),b3MakeVector3(1,0,0),6);
	m_instancingRenderer->drawPoint(b3MakeVector3(0,1,0),b3MakeVector3(0,1,0),6);
	m_instancingRenderer->drawPoint(b3MakeVector3(0,0,1),b3MakeVector3(0,0,1),6);
}

SimpleOpenGL3App::~SimpleOpenGL3App()
{
	delete m_primRenderer ;

	m_window->closeWindow();
	delete m_window;
	delete m_data ;
}


void SimpleOpenGL3App::swapBuffer()
{
	{
		///wireframe line debug rendering
		b3Vector4 lineColor=b3MakeVector4(0,0,1,1);
		b3Vector4 pointColor=b3MakeVector4(1,0,0,1);
		if (m_indices.size())
		{
			m_instancingRenderer->drawLines(&m_vertices[0].x,lineColor,m_vertices.size(),sizeof(b3Vector3),&m_indices[0],m_indices.size(),3);
		}
		m_instancingRenderer->drawPoints(&m_vertices[0].x,pointColor,m_vertices.size(),sizeof(b3Vector3),3);
		m_indices.clear();
		m_vertices.clear();
	}

	m_window->endRendering();
	m_window->startRendering();
}