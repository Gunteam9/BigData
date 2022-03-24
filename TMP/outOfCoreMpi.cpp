#include <string>
#include <iostream>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>

#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPNGWriter.h>
#include <vtkImageData.h>
#include <vtkCamera.h>
#include <vtkContourFilter.h>
#include <vtkProperty.h>
#include <math.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

#include "config.h"
#include "helpers.h"
#include "mpi.h"
#include "globals.h"

#define BIG


const char *location = FICHIER;
int winSize = 500;
int numPasses = 6;
int nbimages = 10;
const char *prefix = "";
int passNum = 0;

int parRank = 0;
int parSize = 1;

using namespace std;

// Function prototypes
vtkRectilinearGrid *ReadGrid(int zStart, int zEnd);

void WriteImage(const char *name, const float *rgba, int width, int height);

vtkRectilinearGrid *ParallelReadGrid() {
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

bool CompositeImageParallel(const float *rgbaIn, float *zBuffer, float *rgbaOut, int imageWidth, int imageHeight) {
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

bool ComposeImageZbuffer(float *rgba_out, float *zbuffer,   int image_width, int image_height);

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &parRank);
    MPI_Comm_size(MPI_COMM_WORLD, &parSize);

    int startexplorevallocal, endexplorevallocal, plage;

    plage = endexploreval - startexploreval;
    if (parRank==0) {
        startexplorevallocal = startexploreval;
        endexplorevallocal = startexplorevallocal+plage/parSize;
    }
    else {
        startexplorevallocal = startexploreval+(plage/parSize)*parRank;
        endexplorevallocal = startexplorevallocal+plage/parSize;
    }

    int interval = 3000; // Should be changed

    for(int i = startexplorevallocal; i < (endexplorevallocal-interval); i += interval){
        vtkContourFilter *cf = vtkContourFilter::New();
        cf->SetNumberOfContours(1);
        int valcont=i;
        cf->SetValue(1,valcont);


        int maxsize=std::max(gridSize,std::max(YgridSize,ZgridSize));
        vtkSmartPointer<vtkTransform> transform =vtkSmartPointer<vtkTransform>::New();
        transform->Scale(gridSize/(float)maxsize,YgridSize/(float)maxsize,ZgridSize/(float)maxsize);
        vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
        transformFilter->SetInputConnection(cf->GetOutputPort());
        transformFilter->SetTransform(transform);

        vtkDataSetMapper *mapper = vtkDataSetMapper::New();
        mapper->SetInputConnection(cf->GetOutputPort());

        vtkLookupTable *lut = vtkLookupTable::New();

        lut->SetHueRange(       0.1, 0.0);
        lut->SetSaturationRange(0.0, 1.0);
        lut->SetValueRange(     1.0, 255.0);
        lut->SetNumberOfColors(100);
        lut->Build();
        mapper->SetLookupTable(lut);
        mapper->SetScalarRange(i,i+interval);

        vtkActor *actor = vtkActor::New();
        actor->SetMapper(mapper);

        vtkRenderer *ren = vtkRenderer::New();
        ren->AddActor(actor);

        vtkCamera* cam = ren->GetActiveCamera();
        cam->SetFocalPoint(0.5, 0.5, 0.5);
        cam->SetPosition(-0., .0, 3.);
        cam->SetViewUp(0., -1.0, 0.0);
        cam->Azimuth(90);

        vtkRenderWindow *renwin = vtkRenderWindow::New();
        renwin->OffScreenRenderingOn();
        renwin->SetSize(winSize, winSize);
        renwin->AddRenderer(ren);

        float *rgba;
        float *auxrgba = new float[4*winSize*winSize];
        float *auxzbuffer = new float[4*winSize*winSize];

        for (int i = 0 ; i < winSize*winSize ; i++)
        {  auxzbuffer[i]=1.0;
            auxrgba[i*4] =  1;
            auxrgba[i*4+1] = 1;
            auxrgba[i*4+2] = 1;
            auxrgba[i*4+3] = 0;
        }


        for(int passCounter=0; passCounter<numPasses;passCounter++){

            int Step=(ZgridSize/numPasses)-1;
            int zStart = passCounter*Step;
            int zEnd = zStart+Step;

            vtkRectilinearGrid *rg = ReadGrid(zStart, zEnd);

            cf->SetInputData(rg);
            rg->Delete();

            cf->Update();
            cf->GetOutput()->GetPointData()->SetActiveScalars("pass_count");

            renwin->Render();

            rgba = renwin->GetRGBAPixelData(0, 0, winSize-1, winSize-1, 1);
            float *zbuffer = renwin->GetZbufferData(0, 0, winSize-1, winSize-1);

            for (int i = 0 ; i < winSize*winSize ; i++)
            {  if (auxzbuffer[i] > zbuffer[i]) {
                    auxzbuffer[i]=zbuffer[i];
                    auxrgba[i*4] = rgba[i*4];
                    auxrgba[i*4+1] = rgba[i*4+1];
                    auxrgba[i*4+2] = rgba[i*4+2];
                    auxrgba[i*4+3] = rgba[i*4+3];
                }
            }


            float *rgba_buffer = new float[4*winSize*winSize];
            bool composer = ComposeImageZbuffer(rgba_buffer, zbuffer,  winSize, winSize);


            free(rgba);
            free(zbuffer);
            free(rgba_buffer);
        }


        char name[128];
        sprintf(name, "final_image_%d_%d.png", i, i+interval);

        WriteImage(name, auxrgba, winSize, winSize);
        free(auxrgba);
        free(auxzbuffer);
    }

    MPI_Finalize();
}



// You should not need to modify these routines.

vtkRectilinearGrid *ReadGrid(int zStart, int zEnd) {
    int i;

    /*   if (zStart < 0 || zEnd < 0 || zStart >= gridSize || zEnd >= gridSize || zStart > zEnd)
     {
     cerr << prefix << "Invalid range: " << zStart << "-" << zEnd << endl;
     return NULL;
     }
     */
    ifstream ifile(location);
    if (ifile.fail()) {
        cerr << prefix << "Unable to open file: " << location << "!!" << endl;
    }

    cerr << prefix << "Reading from " << zStart << " to " << zEnd << endl;

    vtkRectilinearGrid *rg = vtkRectilinearGrid::New();
    vtkFloatArray *X = vtkFloatArray::New();
    X->SetNumberOfTuples(gridSize);
    for (i = 0; i < gridSize; i++)
        X->SetTuple1(i, i / (gridSize - 1.0));
    rg->SetXCoordinates(X);
    X->Delete();
    vtkFloatArray *Y = vtkFloatArray::New();
    Y->SetNumberOfTuples(YgridSize);
    for (i = 0; i < YgridSize; i++)
        Y->SetTuple1(i, i / (YgridSize - 1.0));
    rg->SetYCoordinates(Y);
    Y->Delete();
    vtkFloatArray *Z = vtkFloatArray::New();
    int numSlicesToRead = zEnd - zStart + 1;
    Z->SetNumberOfTuples(numSlicesToRead);
    for (i = zStart; i <= zEnd; i++)
        Z->SetTuple1(i - zStart, i / (ZgridSize - 1.0));
    rg->SetZCoordinates(Z);
    Z->Delete();

    rg->SetDimensions(gridSize, YgridSize, numSlicesToRead);

    unsigned int valuesPerSlice = gridSize * YgridSize;

#if defined(SHORT)
    unsigned int bytesPerSlice = sizeof(short) * valuesPerSlice;

#elif defined(CHAR)
    unsigned int bytesPerSlice   = sizeof(char)*valuesPerSlice;

#elif  defined(FLOAT)
    unsigned int bytesPerSlice   = sizeof(float)*valuesPerSlice;

#else
#error Unsupported choice setting
#endif


#if defined(SMALL)
    unsigned int offset          = (unsigned int)zStart * (unsigned int)bytesPerSlice;
   unsigned int bytesToRead     = bytesPerSlice*numSlicesToRead;
   unsigned int valuesToRead    = valuesPerSlice*numSlicesToRead;
#elif defined(BIG)
    unsigned long long offset = (unsigned long long) zStart * bytesPerSlice;
    unsigned long long bytesToRead = (unsigned long long) bytesPerSlice * numSlicesToRead;
    unsigned int valuesToRead = (unsigned int) valuesPerSlice * numSlicesToRead;
#else
#error Unsupported choice setting
#endif


#if defined(SHORT)
    vtkUnsignedShortArray *scalars = vtkUnsignedShortArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    unsigned short *arr = scalars->GetPointer(0);

#elif defined(CHAR)
    vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    unsigned char *arr = scalars->GetPointer(0);
    
#elif  defined(FLOAT)
    vtkFloatArray *scalars = vtkFloatArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    float *arr = scalars->GetPointer(0);
#else
#error Unsupported choice setting
#endif


    ifile.seekg(offset, ios::beg);
    ifile.read((char *) arr, bytesToRead);
    ifile.close();

    int min = +2147483647;
    int max = 0;


    for (int i = 0; i < valuesToRead; i++) {
        if (min > (scalars->GetPointer(0))[i]) min = (scalars->GetPointer(0))[i];
        if (max < (scalars->GetPointer(0))[i]) max = (scalars->GetPointer(0))[i];

        if (rand() % (valuesToRead / 20) == 0) {
#if defined(SHORT)
            std::cout << (unsigned short) (scalars->GetPointer(0))[i] << " ";

#elif defined(CHAR)
            std::cout<<+(unsigned char)(scalars->GetPointer(0))[i]<<" ";

#elif  defined(FLOAT)
            std::cout<<(float)(scalars->GetPointer(0))[i]<<" ";

#else
#error Unsupported choice setting
#endif


        }
    }


    std::cout << endl << fflush;
    std::cout << "min value read: " << min << endl << fflush;
    std::cout << "max value read: " << max << endl << fflush;


    scalars->SetName("entropy");
    rg->GetPointData()->AddArray(scalars);
    scalars->Delete();

    vtkFloatArray *pr = vtkFloatArray::New();
    pr->SetNumberOfTuples(valuesToRead);
    for (int i = 0; i < valuesToRead; i++)
        pr->SetTuple1(i, passNum);

    pr->SetName("pass_num");
    rg->GetPointData()->AddArray(pr);
    pr->Delete();

    rg->GetPointData()->SetActiveScalars("entropy");

    cerr << prefix << " Done reading" << endl;
    return rg;
}


void WriteImage(const char *name, const float *rgba, int width, int height) {
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