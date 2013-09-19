#ifndef SIMPLE_OPENGL3_APP_H
#define SIMPLE_OPENGL3_APP_H

#include "OpenGLWindow/GLInstancingRenderer.h"
#include "OpenGLWindow/GLPrimitiveRenderer.h"
#include "OpenGLWindow/b3gWindowInterface.h"
#include "../btgui/Bullet3Common/b3AlignedObjectArray.h"
class b3Vector3;

struct SimpleOpenGL3App
{
	struct SimpleInternalData* m_data;

	class b3gWindowInterface*	m_window;
	class GLPrimitiveRenderer*	m_primRenderer;
	class GLInstancingRenderer* m_instancingRenderer;
	b3AlignedObjectArray<unsigned int>	m_indices;
	b3AlignedObjectArray<b3Vector3>		m_vertices;


	void drawArc(const b3Vector3& center, const b3Vector3& normal, const b3Vector3& axis, b3Scalar radiusA, b3Scalar radiusB, b3Scalar minAngle, b3Scalar maxAngle,  bool drawSect, b3Scalar stepDegrees = b3Scalar(10.f));
	
public:
	SimpleOpenGL3App(const char* title, int width,int height);
	virtual ~SimpleOpenGL3App();
	
	
	void	swapBuffer();

	void	drawText( const char* txt, int posX, int posY);

	void	drawGrid(int gridSize=10);
	///drawing of some basic wireframe shapes
	void	drawLine(const float* from,const float* to);
	void	drawCylinder(float radius, float halfHeight, int upAxis, const float* pos, const float* ornQuat);
	void	drawCone(float radius, float height, int upAxis, const float* pos, const float* ornQuat);
	void	drawSpherePatch(const b3Vector3& center, const b3Vector3& up, const b3Vector3& axis, b3Scalar radius, 
							b3Scalar minTh, b3Scalar maxTh, b3Scalar minPs, b3Scalar maxPs, b3Scalar stepDegrees = b3Scalar(10.f),bool drawCenter = true);
	void	drawCapsule(float radius, float halfHeight, int upAxis,  const float* pos, const float* ornQuat);
	void	drawSphere(float radius,  const float* pos, const float* ornQuat);
	void	drawBox(const float* halfDimensions, const float* position, const float* orientation);
};

#endif //SIMPLE_OPENGL3_APP_H
