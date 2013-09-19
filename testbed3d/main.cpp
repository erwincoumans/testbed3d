
#include "Bullet3Common/b3Vector3.h"
#include "SimpleOpenGL3App.h"
#include "Bullet3Common/b3Transform.h"







int main(int argc, char* argv[])
{
	
	//float dt = 1./120.f;
	float dt = 1./30.f;
	
	SimpleOpenGL3App* app = new SimpleOpenGL3App("title",1024,768);
	app->m_instancingRenderer;

	app->m_instancingRenderer->setCameraDistance(13);
	app->m_instancingRenderer->setCameraPitch(0);
	app->m_instancingRenderer->setCameraTargetPosition(b3MakeVector3(0,0,0));

	GLint err = glGetError();
  b3Assert(err==GL_NO_ERROR);




	do
	{
		
		GLint err = glGetError();
		b3Assert(err==GL_NO_ERROR);
		app->m_instancingRenderer->init();
		app->m_instancingRenderer->updateCamera();
		

		app->drawGrid();
		app->drawText("Text test",10,10);
		
		b3Vector3 pos=b3MakeVector3(0,0,0);
		b3Quaternion orn(0,0,0,1);
		b3Vector3 size=b3MakeVector3(1,1,1,1);


		app->drawBox(size, pos,orn);
		b3Transform tr; tr.setIdentity();
		tr.setOrigin(b3MakeVector3(2,0,0));

		app->drawCapsule(1,2,1,tr.getOrigin(),tr.getRotation());
		tr.setOrigin(b3MakeVector3(4,0,0));
		app->drawSphere(1,tr.getOrigin(),tr.getRotation());
		tr.setOrigin(b3MakeVector3(6,0,0));
		app->drawCone(1,1,1,tr.getOrigin(),tr.getRotation());
		tr.setOrigin(b3MakeVector3(8,0,0));
		app->drawCylinder(1,1,2,tr.getOrigin(),tr.getRotation());
		

				
	
		
		
		app->swapBuffer();
	} while (!app->m_window->requestedExit());

	
	delete app;
	return 0;
}
