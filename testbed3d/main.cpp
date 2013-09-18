
#include "Bullet3Common/b3Vector3.h"
#include "SimpleOpenGL3App.h"
#include "Bullet3Common/b3Transform.h"



b3AlignedObjectArray<unsigned int> indices;
b3AlignedObjectArray<b3Vector3> vertices;
class GLInstancingRenderer* gRenderer=0;

void	drawLine(const float* from,const float* to)
{
	
	int numV = vertices.size();
	vertices.push_back(b3MakeVector3(from[0],from[1],from[2]));
	vertices.push_back(b3MakeVector3(to[0],to[1],to[2]));
	indices.push_back(numV);
	indices.push_back(numV+1);
}




void drawArc(const b3Vector3& center, const b3Vector3& normal, const b3Vector3& axis, b3Scalar radiusA, b3Scalar radiusB, b3Scalar minAngle, b3Scalar maxAngle, 
			 bool drawSect, b3Scalar stepDegrees = b3Scalar(10.f))
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
	

void drawCylinder(float radius, float halfHeight, int upAxis, const float* pos, const float* ornQuat)
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

void drawCone(float radius, float height, int upAxis, const float* pos, const float* ornQuat)
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

void drawSpherePatch(const b3Vector3& center, const b3Vector3& up, const b3Vector3& axis, b3Scalar radius, 
		b3Scalar minTh, b3Scalar maxTh, b3Scalar minPs, b3Scalar maxPs, b3Scalar stepDegrees = b3Scalar(10.f),bool drawCenter = true)
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

void drawCapsule(float radius, float halfHeight, int upAxis,  const float* pos, const float* ornQuat)
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

 void	drawSphere(float radius,  const float* pos, const float* ornQuat)
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


void drawBox(const float* halfDimensions, const float* position, const float* orientation)
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

	int startIndex=vertices.size();

	vertices.push_back(org - dx - dy - dz);
	vertices.push_back(org + dx - dy - dz);
	vertices.push_back(org + dx - dy + dz);
	vertices.push_back(org - dx - dy + dz);
	vertices.push_back(org - dx + dy - dz);
	vertices.push_back(org + dx + dy - dz);
	vertices.push_back(org + dx + dy + dz);
	vertices.push_back(org - dx + dy + dz);
	
	float lineColor[4]={0.,0.,.1,1};

	float lineSize = 2;
	
	
	indices.push_back(startIndex+0); indices.push_back(startIndex+1);
	indices.push_back(startIndex+1); indices.push_back(startIndex+2);
	indices.push_back(startIndex+2); indices.push_back(startIndex+3);
	indices.push_back(startIndex+3); indices.push_back(startIndex+0);

	indices.push_back(startIndex+4); indices.push_back(startIndex+5);
	indices.push_back(startIndex+5); indices.push_back(startIndex+6);
	indices.push_back(startIndex+6); indices.push_back(startIndex+7);
	indices.push_back(startIndex+7); indices.push_back(startIndex+4);

	indices.push_back(startIndex+0); indices.push_back(startIndex+4);
	indices.push_back(startIndex+1); indices.push_back(startIndex+5);
	indices.push_back(startIndex+2); indices.push_back(startIndex+6);
	indices.push_back(startIndex+3); indices.push_back(startIndex+7);

	

}


int main(int argc, char* argv[])
{
	
	//float dt = 1./120.f;
	float dt = 1./30.f;
	
	SimpleOpenGL3App* app = new SimpleOpenGL3App("title",1024,768);
	gRenderer = app->m_instancingRenderer;

	app->m_instancingRenderer->setCameraDistance(13);
	app->m_instancingRenderer->setCameraPitch(0);
	app->m_instancingRenderer->setCameraTargetPosition(b3MakeVector3(0,0,0));

	GLint err = glGetError();
  b3Assert(err==GL_NO_ERROR);




	do
	{
		vertices.clear();
		indices.clear();

		GLint err = glGetError();
		b3Assert(err==GL_NO_ERROR);
		app->m_instancingRenderer->init();
		app->m_instancingRenderer->updateCamera();
		

		app->drawGrid();
		app->drawText("Text test",10,10);
		
		b3Vector3 pos=b3MakeVector3(0,0,0);
		b3Quaternion orn(0,0,0,1);
		b3Vector3 size=b3MakeVector3(1,1,1,1);


		drawBox(size, pos,orn);
		b3Transform tr; tr.setIdentity();
		tr.setOrigin(b3MakeVector3(2,0,0));

		drawCapsule(1,2,1,tr.getOrigin(),tr.getRotation());
		tr.setOrigin(b3MakeVector3(4,0,0));
		drawSphere(1,tr.getOrigin(),tr.getRotation());
		tr.setOrigin(b3MakeVector3(6,0,0));
		drawCone(1,1,1,tr.getOrigin(),tr.getRotation());
		tr.setOrigin(b3MakeVector3(8,0,0));
		drawCylinder(1,1,2,tr.getOrigin(),tr.getRotation());
		

				
		b3Vector4 lineColor=b3MakeVector4(0,0,1,1);
		b3Vector4 pointColor=b3MakeVector4(1,0,0,1);
		
		if (indices.size())
		{
			gRenderer->drawLines(&vertices[0].x,lineColor,vertices.size(),sizeof(b3Vector3),&indices[0],indices.size(),3);
		}
		//gRenderer->drawPoints(&vertices[0].x,pointColor,vertices.size(),sizeof(b3Vector3),6);
		
		indices.clear();
		vertices.clear();
		app->swapBuffer();
	} while (!app->m_window->requestedExit());

	
	delete app;
	return 0;
}
