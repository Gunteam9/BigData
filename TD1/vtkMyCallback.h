//
// Created by o2196290 on 02/03/2022.
//

#ifndef VTKTP_VTKMYCALLBACK_H
#define VTKTP_VTKMYCALLBACK_H

#include <vtkCommand.h>
#include <vtkProperty.h>
#include <vtkBoxWidget.h>
#include <vtkTransform.h>
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


class vtkMyCallback : public vtkCommand
{
    public:
    static vtkMyCallback* New()
    {
        return new vtkMyCallback;
    }

        virtual void Execute(vtkObject *caller, unsigned long, void*);
};


#endif //VTKTP_VTKMYCALLBACK_H
