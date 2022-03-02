#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProgrammableFilter.h>
#include <vtkCallbackCommand.h>
#include <vtkConeSource.h>
#include <vtkProperty.h>
#include "vtkCamera.h"
#include "unistd.h"
#include "vtkMyCallback.h"

using namespace std;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "EndlessLoop"

// Globals
unsigned int counter = 0;
const unsigned int coneCount = 5;

void TimerCallbackFunction ( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData) )
{
    cout << "timer callback" << endl;
    vtkSmartPointer<vtkRenderer> renderer = static_cast<vtkRenderer*>(clientData);
//    renderer->GetActiveCamera()->Azimuth( 10 );

    vtkActorCollection* actorCollection = renderer->GetActors();
    vtkCollectionSimpleIterator actorIte;
    vtkActor* currentActor;
    actorCollection->InitTraversal(actorIte);
    while (currentActor = actorCollection->GetNextActor(actorIte)) {
        currentActor->GetProperty()->SetColor((float) rand() / RAND_MAX, (float) rand() / RAND_MAX, (float) rand() / RAND_MAX);
    }

    vtkRenderWindowInteractor* iren = dynamic_cast<vtkRenderWindowInteractor*>(caller);
    iren->GetRenderWindow()->Render();

    iren->Render();

    counter++;
}

int main(int, char *[])
{

    srand(time(NULL));

    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetCenter(0.0, 0.0, 0.0);
    sphereSource->SetRadius(5.0);
    sphereSource->Update();

    vtkConeSource** cones = new vtkConeSource*[coneCount];
    for (int i = 0; i < coneCount; ++i) {
        cones[i] = vtkConeSource::New();
        cones[i]->SetHeight(2);
        cones[i]->SetRadius(0.5);
        cones[i]->SetResolution(50);
        cones[i]->Update();
    }

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(cones[0]->GetOutputPort());
    for (int i = 1; i < coneCount; ++i)
        mapper->AddInputConnection(cones[i]->GetOutputPort());

    // Create an actor
    vtkSmartPointer<vtkActor>* actors = new vtkSmartPointer<vtkActor>[coneCount];
    for (int i = 0; i < coneCount; ++i) {
        actors[i] = vtkSmartPointer<vtkActor>::New();
        actors[i]->SetPosition(rand() % 10,rand() % 10,rand() % 10);
        actors[i]->SetMapper(mapper);

        vtkProperty* property = vtkProperty::New();
        actors[i]->SetProperty(property);
    }

    // Setup renderer, render window, and interactor
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

    for (int i = 0; i < coneCount; ++i)
        renderer->AddActor(actors[i]);
    renderer->SetBackground(0, 0, 0);

    vtkSmartPointer<vtkRenderWindow> renderWindow =vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> iren =vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renderWindow);
    iren->Initialize();
    iren->CreateRepeatingTimer(100);

    vtkSmartPointer<vtkCallbackCommand> timerCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    timerCallback->SetCallback(TimerCallbackFunction);
    timerCallback->SetClientData(renderer);

    iren->AddObserver ( vtkCommand::TimerEvent, timerCallback );
//    iren->Start();

    for (int i = 0; i < coneCount; ++i) {
        vtkBoxWidget *boxWidget = vtkBoxWidget::New();
        boxWidget->SetInteractor(iren);
        boxWidget->SetPlaceFactor(1.25);
        boxWidget->SetProp3D(actors[i]);
        boxWidget->PlaceWidget();
        vtkMyCallback* callback = vtkMyCallback::New();
        boxWidget->AddObserver(vtkCommand::InteractionEvent, callback);
        boxWidget->On();
    }

    iren->Initialize();
    iren->Start();

    // Start a timer 10 seconds to keep visible the rendering Window
//    usleep(10 * 1000000);

    return EXIT_SUCCESS;
}

#pragma clang diagnostic pop