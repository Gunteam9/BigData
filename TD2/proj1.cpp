
#include <vtkPolyData.h>
#include <vtkPLYReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>

#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>


#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleTerrain.h>
#include <vtkInteractorStyleJoystickCamera.h>
#include <vtkInteractorStyleFlight.h>


#include <vtkElevationFilter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

#include "vtkCamera.h"
#include "vtkProperty.h"

#include "vtkSphereSource.h"
#include "vtkConeSource.h"

#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkContourFilter.h>
#include <vtkCutter.h>
#include <vtkPlane.h>


int main(int argc, char **argv)
{
    vtkDataSetReader *reader = vtkDataSetReader::New();
    reader->SetFileName("../noise.vtk");

    vtkContourFilter *cf = vtkContourFilter::New();
    cf->SetNumberOfContours(2);
    cf->SetValue(0, 2.4);
    cf->SetValue(1, 4.0);
    cf->SetInputConnection(reader->GetOutputPort());


    vtkCutter* cutter = vtkCutter::New();
    vtkPlane* plane = vtkPlane::New();
    plane->SetNormal(0,1,1);
    plane->SetOrigin(0,0,0);
    cutter->SetCutFunction(plane);
    cutter->SetInputConnection(reader->GetOutputPort());

    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
//    mapper->SetInputConnection(cutter->GetOutputPort());
    mapper->SetInputConnection(cf->GetOutputPort());

    vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);
    
    vtkRenderer *ren = vtkRenderer::New();
    ren->SetViewport(0.5, 0, 1, 1);
    ren->AddActor(actor);

    vtkRenderer *ren2 = vtkRenderer::New();
    ren2->SetViewport(0, 0, 0.5, 1);



    vtkLookupTable *lut = vtkLookupTable::New();
//    mapper->SetLookupTable(lut);
//    mapper->SetScalarRange(1,6);
//    lut->Build();

    for (int i = 0; i < 256; ++i) {
        double vals[4] = {(255 - i) / 255.0, 0, i / 255.0, 1};
        lut->SetTableValue(i, vals);
    }
    lut->SetNumberOfColors(256);
    mapper->SetLookupTable(lut);
    mapper->SetScalarRange(2.4,6);
    lut->Build();

    vtkDataSetMapper *mapper2 = vtkDataSetMapper::New();
    vtkActor *actor2 = vtkActor::New();
    actor2->SetMapper(mapper2);
    mapper2->SetInputConnection(cutter->GetOutputPort());
    mapper2->SetScalarRange(1,6);
    mapper2->SetLookupTable(lut);
    ren2->AddActor(actor2);

    vtkRenderWindow *renwin = vtkRenderWindow::New();
    renwin->SetSize(768, 768);
    renwin->AddRenderer(ren);
    renwin->AddRenderer(ren2);

    int nbpas=200;
    for (int i = 0 ; i < nbpas ; i++)
    {
        cf->SetValue(i, (((float)i/nbpas) * 5) + 1);
        renwin->Render();
    }

//    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
//    iren->SetRenderWindow(renwin);
//    renwin->Render();
//    iren->Start();
}

