//
// Created by o2196290 on 02/03/2022.
//

#include "vtkMyCallback.h"

void vtkMyCallback::Execute(vtkObject *caller, unsigned long, void *) {
    vtkTransform *t = vtkTransform::New();
    vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
    widget->GetTransform(t);
    widget->GetProp3D()->SetUserTransform(t);
    t->Delete();
}
