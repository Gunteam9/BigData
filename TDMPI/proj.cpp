
#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkContourFilter.h>
#include <vtkPointData.h>
#include <vtkGraphicsFactory.h>
//#include <vtkImagingFactory.h>
#include <vtkPNGWriter.h>
#include <vtkImageData.h>
#include <vtkCamera.h>
#include <unistd.h>
#include <chrono>
#include "config.h"
#include "helpers.h"

#define STR_VALUE(arg) #arg
#define FUNCTION_NAME(name) STR_VALUE(name)
#define TAILLE_NAME FUNCTION_NAME(TAILLE)

//ICI PLACER LA TAILLE DU FICHIER 512 ou 1024
#define TAILLE 512
#define FICHIER MY_MESHES_PATH "/sn_resamp" TAILLE_NAME

int winSize = 768;
const char *prefix = "";;
int gridSize = TAILLE;
const char *location = FICHIER;

int parRank = 0;
int parSize = 1;

using std::cerr;
using std::endl;

// Function prototypes
vtkRectilinearGrid *ReadGrid(int zStart, int zEnd);

void WriteImage(const char *name, const float *rgba, int width, int height);


vtkRectilinearGrid* ParallelReadGrid() {
    int zCount = parRank < gridSize % parSize
                 ? (gridSize / parSize) + 1
                 : gridSize / parSize;
    int zStart = parRank < gridSize % parSize
                 ? ((gridSize / parSize) + 1) * parRank
                 : ((gridSize / parSize) + 1) * (gridSize % parSize) +
                   (gridSize / parSize) * (parRank - (gridSize % parSize));
    int zEnd = zStart + zCount;

    if (parRank == parSize - 1)
        zEnd -= 1;

    return ReadGrid(zStart, zEnd);
}

bool CompositeImage(const float *rgbaIn, float *zBuffer, float *rgbaOut, int imageWidth, int imageHeight) {
    int nPixels = imageWidth * imageHeight;
    auto *zBufferMin = new float[nPixels];

    MPI_Allreduce(zBuffer, zBufferMin, nPixels, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

    auto *rgbaTmp = new float[4 * nPixels];
    for (int i = 0; i < nPixels; ++i) {
        if (zBuffer[i] <= zBufferMin[i]) {
            rgbaTmp[i * 4 + 0] = rgbaIn[i * 4 + 0];
            rgbaTmp[i * 4 + 1] = rgbaIn[i * 4 + 1];
            rgbaTmp[i * 4 + 2] = rgbaIn[i * 4 + 2];
            rgbaTmp[i * 4 + 3] = rgbaIn[i * 4 + 3];
        }
    }

    MPI_Reduce(rgbaTmp, rgbaOut, 4 * nPixels, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

    delete[] zBufferMin;
    delete[] rgbaTmp;

    return true;
}

int main(int argc, char *argv[]) {
// MPI setup
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &parRank);
    MPI_Comm_size(MPI_COMM_WORLD, &parSize);

// I use the variable "prefix" to get the print statements right.
    std::string p(std::to_string(parRank));
    int t1;
    t1 = timer->StartTimer();


    prefix = p.c_str();

    GetMemorySize((p + ":initialization").c_str());


// Read the data.
    vtkRectilinearGrid *rg = ParallelReadGrid();
    GetMemorySize((p + ":After Read").c_str());

// Contour the data.
    vtkContourFilter *cf = vtkContourFilter::New();
    cf->SetNumberOfContours(1);
    cf->SetValue(0, 3.4);
    cf->SetInputData(rg);

// Force an update and set the parallel rank as the active scalars.
    cf->Update();
    cf->GetOutput()->GetPointData()->SetActiveScalars("par_rank");

    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
    mapper->SetInputConnection(cf->GetOutputPort());

    vtkLookupTable *lut = vtkLookupTable::New();
    mapper->SetLookupTable(lut);
    mapper->SetScalarRange(0, parSize - 1.0);
    lut->Build();

    vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);

    vtkRenderer *ren = vtkRenderer::New();
    ren->AddActor(actor);

    vtkCamera *cam = ren->GetActiveCamera();
    cam->SetFocalPoint(0.5, 0.5, 0.5);
    cam->SetPosition(1.5, 1.5, 1.5);

    vtkRenderWindow *renwin = vtkRenderWindow::New();
// THIS DOESN'T WORK WITHOUT MESA
    renwin->OffScreenRenderingOn();
    renwin->SetSize(winSize, winSize);
    renwin->AddRenderer(ren);


//WE'RE JUST RENDERING A SINGLE IMAGE ... NO INTERACTOR
/*   vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
   iren->SetRenderWindow(renwin);
   renwin->Render();
   iren->Start();
 */

    float *rgba, *zbuffer;
    bool staggerGLRenders = false;
//   bool staggerGLRenders = true;

    if (staggerGLRenders) {
        for (int i = 0; i < parSize; i++) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (i == parRank) {
                renwin->Render();
                sleep(1);
                rgba = renwin->GetRGBAPixelData(0, 0, (winSize - 1), (winSize - 1), 1);
                zbuffer = renwin->GetZbufferData(0, 0, (winSize - 1), (winSize - 1));
            }
        }
    } else {
        renwin->Render();
        rgba = renwin->GetRGBAPixelData(0, 0, (winSize - 1), (winSize - 1), 1);
        zbuffer = renwin->GetZbufferData(0, 0, (winSize - 1), (winSize - 1));
    }

    float *new_rgba = new float[4 * winSize * winSize];
    bool didComposite = CompositeImage(rgba, zbuffer, new_rgba, winSize, winSize);
    if (didComposite) {
        if (parRank == 0) {
            WriteImage("final_image.png", new_rgba, winSize, winSize);
        }

        char name[128];
        sprintf(name, "img%d.png", parRank);
        WriteImage(name, rgba, winSize, winSize);

    }

    renwin->Delete();
    rg->Delete();
    GetMemorySize((p + ":End").c_str());
    timer->StopTimer(t1, (p + ":Time").c_str());
    MPI_Finalize();
}


// You should not need to modify these routines.


// You should not need to modify these routines.
vtkRectilinearGrid *
ReadGrid(int zStart, int zEnd) {
    int i;
    if (zStart < 0 || zEnd < 0 || zStart >= gridSize || zEnd >= gridSize || zStart > zEnd) {
        cerr << prefix << "Invalid range: " << zStart << "-" << zEnd << endl;
        return NULL;
    }

    ifstream ifile(location);
    if (ifile.fail()) {
        cerr << prefix << "Unable to open file: " << location << "!!" << endl;
    }

    cerr << prefix << "Reading " << location << " from " << zStart << " to " << zEnd << endl;

    vtkRectilinearGrid *rg = vtkRectilinearGrid::New();
    vtkFloatArray *X = vtkFloatArray::New();
    X->SetNumberOfTuples(gridSize);
    for (i = 0; i < gridSize; i++)
        X->SetTuple1(i, i / (gridSize - 1.0));
    rg->SetXCoordinates(X);
    X->Delete();
    vtkFloatArray *Y = vtkFloatArray::New();
    Y->SetNumberOfTuples(gridSize);
    for (i = 0; i < gridSize; i++)
        Y->SetTuple1(i, i / (gridSize - 1.0));
    rg->SetYCoordinates(Y);
    Y->Delete();
    vtkFloatArray *Z = vtkFloatArray::New();
    int numSlicesToRead = zEnd - zStart + 1;
    Z->SetNumberOfTuples(numSlicesToRead);
    for (i = zStart; i <= zEnd; i++)
        Z->SetTuple1(i - zStart, i / (gridSize - 1.0));
    rg->SetZCoordinates(Z);
    Z->Delete();

    rg->SetDimensions(gridSize, gridSize, numSlicesToRead);


    int valuesPerSlice = gridSize * gridSize;
    int bytesPerSlice = 4 * valuesPerSlice;


#if TAILLE == 512
    unsigned int offset = (unsigned int) zStart * (unsigned int) bytesPerSlice;
    unsigned int bytesToRead = bytesPerSlice * numSlicesToRead;
    unsigned int valuesToRead = valuesPerSlice * numSlicesToRead;
#elif TAILLE == 1024
    unsigned long long offset = (unsigned long long) zStart * bytesPerSlice;
    unsigned long long bytesToRead = (unsigned long long) bytesPerSlice * numSlicesToRead;
    unsigned int valuesToRead = (unsigned int) valuesPerSlice * numSlicesToRead;
#else
#error Unsupported choice setting
#endif


    vtkFloatArray *scalars = vtkFloatArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    float *arr = scalars->GetPointer(0);
    ifile.seekg(offset, ios::beg);
    ifile.read((char *) arr, bytesToRead);
    ifile.close();

    scalars->SetName("entropy");
    rg->GetPointData()->AddArray(scalars);
    scalars->Delete();


    vtkFloatArray *pr = vtkFloatArray::New();
    pr->SetNumberOfTuples(valuesToRead);
    for (int i = 0; i < valuesToRead; i++)
        pr->SetTuple1(i, parRank);

    pr->SetName("par_rank");
    rg->GetPointData()->AddArray(pr);
    pr->Delete();

    rg->GetPointData()->SetActiveScalars("entropy");

    cerr << prefix << " Done reading" << endl;
    return rg;
}

void
WriteImage(const char *name, const float *rgba, int width, int height) {
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
#if VTK_MAJOR_VERSION <= 5
    img->SetNumberOfScalarComponents(3);
    img->SetScalarTypeToUnsignedChar();
#else
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
#endif

    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++) {
            unsigned char *ptr = (unsigned char *) img->GetScalarPointer(i, j, 0);
            int idx = j * width + i;
            ptr[0] = (unsigned char) (255 * rgba[4 * idx + 0]);
            ptr[1] = (unsigned char) (255 * rgba[4 * idx + 1]);
            ptr[2] = (unsigned char) (255 * rgba[4 * idx + 2]);
        }


    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(name);
    writer->Write();

    img->Delete();
    writer->Delete();
}

bool ComposeImageZbuffer(float *rgba_out, float *zbuffer,   int image_width, int image_height)
{
    int npixels = image_width*image_height;

    float min=1;
    float max=0;
    for (int i = 0 ; i < npixels ; i++){
        if (zbuffer[i]<min) min=zbuffer[i];
        if (zbuffer[i]>max) max=zbuffer[i];

    }
    std::cout<<"min:"<<min;
    std::cout<<"max:"<<max<<"  ";

    float coef=1.0/((max-min));

    std::cout<<"coef:"<<coef<<"  ";

    for (int i = 0 ; i < npixels ; i++){

        rgba_out[i*4] = (zbuffer[i]==1.0?0:1-coef*(zbuffer[i]-min));
        rgba_out[i*4+1] = (zbuffer[i]==1.0?0:1-coef*(zbuffer[i]-min));
        rgba_out[i*4+2] = (zbuffer[i]==1.0?0:1-coef*(zbuffer[i]-min));
        rgba_out[i*4+3] = 0.0;
    }


    return false;
}
