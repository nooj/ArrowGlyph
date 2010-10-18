/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkArrowGlyphFilter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkArrowGlyphFilter.h"

#include "vtkCellData.h"
#include "vtkCell.h"
#include "vtkDataSet.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"
#include "vtkArrowSource.h"
#include "vtkMaskPoints.h"
#include "vtkMultiProcessController.h"
#include "vtkProcessModule.h"
#include "vtkTransform.h"

vtkStandardNewMacro(vtkArrowGlyphFilter);
vtkCxxSetObjectMacro(vtkArrowGlyphFilter, ArrowSourceObject, vtkArrowSource);

//----------------------------------------------------------------------------
// Construct object with scaling on, scaling mode is by scalar value,
// scale factor = 1.0, the range is (0,1), orient geometry is on, and
// orientation is by vector. Clamping and indexing are turned off. No
// initial sources are defined.
vtkArrowGlyphFilter::vtkArrowGlyphFilter()
{
  this->Orient        = 1;
  this->OrientArray   = NULL;
  //
  this->LengthScaling = 1;
  this->LengthFactor  = 1.0;
  this->LengthArray   = NULL;
  //
  this->RadiusScaling = 1;
  this->RadiusFactor  = 1.0;
  this->RadiusArray   = NULL;
  //
  this->MaskPoints = vtkMaskPoints::New();
  this->RandomMode = this->MaskPoints->GetRandomMode();
  this->MaximumNumberOfPoints = 5000;
//  this->NumberOfProcesses = vtkMultiProcessController::GetGlobalController() ?
//    vtkMultiProcessController::GetGlobalController()->GetNumberOfProcesses() : 1;
  this->UseMaskPoints = 1;
  //
  this->GeneratePointIds = 0;
  this->PointIdsName = NULL;
  this->SetPointIdsName("InputPointIds");
  this->SetNumberOfInputPorts(1);
  //
  this->ArrowSourceObject = NULL; // vtkSmartPointer<vtkArrowSource>::New();
}

//----------------------------------------------------------------------------
vtkArrowGlyphFilter::~vtkArrowGlyphFilter()
{
  if (this->PointIdsName)
    {
    delete []PointIdsName;
    }
  if (this->OrientArray)
    {
    delete []OrientArray;
    }
  if (this->LengthArray)
    {
    delete []LengthArray;
    }
  if (this->RadiusArray)
    {
    delete []RadiusArray;
    }
  if(this->MaskPoints)
    {
    this->MaskPoints->Delete();
    }
  this->SetArrowSourceObject(NULL);
}

//----------------------------------------------------------------------------
unsigned long vtkArrowGlyphFilter::GetMTime()
{
  unsigned long mTime=this->Superclass::GetMTime();
  unsigned long time;
  if ( this->ArrowSourceObject != NULL )
    {
    time = this->ArrowSourceObject ->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }
  return mTime;
}

//-----------------------------------------------------------------------------
void vtkArrowGlyphFilter::SetRandomMode(int mode)
{
  if (mode==this->MaskPoints->GetRandomMode())
    {
    // no change
    return;
    }
  // Store random mode to so that we don't have to call
  // MaskPoints->GetRandomMode() in tight loop.
  this->MaskPoints->SetRandomMode(mode);
  this->RandomMode = mode;
  this->Modified();
}

//-----------------------------------------------------------------------------
int vtkArrowGlyphFilter::GetRandomMode()
{
  return this->MaskPoints->GetRandomMode();
}

//-----------------------------------------------------------------------------
void vtkArrowGlyphFilter::SetUseMaskPoints(int useMaskPoints)
{
  if (useMaskPoints==this->UseMaskPoints)
    {
    return;
    }
  this->UseMaskPoints=useMaskPoints;
  this->Modified();
}

//-----------------------------------------------------------------------------
vtkIdType vtkArrowGlyphFilter::GatherTotalNumberOfPoints(vtkIdType localNumPts)
{
  // Although this is not perfectly process invariant, it is better
  // than we had before (divide by number of processes).
  vtkIdType totalNumPts = localNumPts;
  vtkMultiProcessController *controller = 
    vtkMultiProcessController::GetGlobalController();
  if (controller)
    {
    vtkIdType tmp;
    // This could be done much easier with MPI specific calls.
    if (controller->GetLocalProcessId() == 0)
      {
      int i;
      // Sum points on all processes.
      for (i = 1; i < controller->GetNumberOfProcesses(); ++i)
        {
        controller->Receive(&tmp, 1, i, vtkProcessModule::GlyphNPointsGather);
        totalNumPts += tmp;
        }
      // Send results back to all processes.
      for (i = 1; i < controller->GetNumberOfProcesses(); ++i)
        {
        controller->Send(&totalNumPts, 1, 
                         i, vtkProcessModule::GlyphNPointsScatter);
        }
      }
    else
      {
      controller->Send(&localNumPts, 1, 
                       0, vtkProcessModule::GlyphNPointsGather);
      controller->Receive(&totalNumPts, 1, 
                          0, vtkProcessModule::GlyphNPointsScatter);
      }
    }

  return totalNumPts;
}

//----------------------------------------------------------------------------
int vtkArrowGlyphFilter::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataObject   *input = inInfo->Get(vtkDataObject::DATA_OBJECT());
  vtkDataSet    *dsInput = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!dsInput) {
    if (input) {
      vtkErrorMacro("This filter cannot process input of type: " << input->GetClassName());
    }
    return 0;
  }

  // Glyph a subset.
  vtkIdType maxNumPts = this->MaximumNumberOfPoints;
  vtkIdType numPts = dsInput->GetNumberOfPoints();
  vtkIdType totalNumPts = this->GatherTotalNumberOfPoints(numPts);

  // What fraction of the points will this processes get allocated?
  maxNumPts = (vtkIdType)((double)(maxNumPts)*(double)(numPts)/(double)(totalNumPts));

  maxNumPts = (maxNumPts < 1) ? 1 : maxNumPts;

  vtkInformationVector* inputVs[2];

  vtkInformationVector* inputV = inputVector[0];
  inputVs[0] = vtkInformationVector::New();
  inputVs[0]->SetNumberOfInformationObjects(1);
  vtkInformation* newInInfo = vtkInformation::New();
  newInInfo->Copy(inputV->GetInformationObject(0));
  inputVs[0]->SetInformationObject(0, newInInfo);
  newInInfo->Delete();
  inputVs[1] = inputVector[1];

  int retVal = this->MaskAndExecute(numPts, maxNumPts, dsInput,
    request, inputVs, outputVector);

  inputVs[0]->Delete();
  return retVal;
}

//----------------------------------------------------------------------------
int vtkArrowGlyphFilter::MaskAndExecute(vtkIdType numPts, vtkIdType maxNumPts,
                                     vtkDataSet* input,
                                     vtkInformation* request,
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector)

{
  //
  // shallow copy input so that internal pipeline doesn't trash information
  // pass input into maskfilter and update
  //
  vtkDataSet* inputCopy = input->NewInstance();
  inputCopy->ShallowCopy(input);
  this->MaskPoints->SetInput(inputCopy);
  inputCopy->Delete();

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  this->MaskPoints->SetMaximumNumberOfPoints(maxNumPts);
  this->MaskPoints->SetOnRatio(numPts / maxNumPts);

  vtkInformation *maskPointsInfo =
    this->MaskPoints->GetExecutive()->GetOutputInformation(0);
  maskPointsInfo->Set(
    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()));
  maskPointsInfo->Set(
    vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()));
  maskPointsInfo->Set(
    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
    outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()));
  this->MaskPoints->Update();
  
  // How many points will we be glyphing (in this process)
  vtkPoints *maskedpoints = this->MaskPoints->GetOutput()->GetPoints();
  vtkIdType Nm = maskedpoints->GetNumberOfPoints();


  //
  // Now we insert the new code specially for our arrow filter
  //

  //
  // get the input and output
  //
  vtkDataSet *minput = this->MaskPoints->GetOutput();

  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  //
  // Make sure our arrow is any good
  //
  if (!this->ArrowSourceObject) {
    vtkArrowSource *arrow = vtkArrowSource::New();
    this->SetArrowSourceObject(arrow);
    arrow->Delete();
  }
  // and get useful information from it
  this->ArrowSourceObject->Update();
  vtkPolyData *arrow = this->ArrowSourceObject->GetOutput();
  vtkPoints *apoints = arrow->GetPoints();
  vtkIdType Na = apoints->GetNumberOfPoints();

  // 
  // Find the arrays to be used for Length/Radius/etc
  // if not present, we will use default values based on particle size
  //
  vtkDataArray *orientdata = this->OrientArray ? minput->GetPointData()->GetArray(this->OrientArray) : NULL;
  vtkDataArray *lengthdata = this->LengthArray ? minput->GetPointData()->GetArray(this->LengthArray) : NULL;
  vtkDataArray *radiusdata = this->RadiusArray ? minput->GetPointData()->GetArray(this->RadiusArray) : NULL;
  bool radiusMagnitude = false;
  bool lengthMagnitude = false;
  if (radiusdata && radiusdata->GetNumberOfComponents()==3) {
    radiusMagnitude = true;
  }
  if (lengthdata && lengthdata->GetNumberOfComponents()==3) {
    lengthMagnitude = true;
  }

  // we know the output will require NumPoints in Arrow * NumPoints in MaskPoints
  // so we can pre-allocate the output space.
  vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
  newPoints->SetNumberOfPoints(Na*Nm);

  vtkTransform *trans = vtkTransform::New();
  

  // The variables we use to control each individual glyph
  double radius = 1.0;
  double length = 1.0;
  double *orientvector = NULL;

  // 
  // Loop over all out points and do the actual glyphing
  //
  for (vtkIdType i=0; i<Nm; i++) {
    
    // update progress bar
    if ( ! (i % 10000) ) {
      this->UpdateProgress(static_cast<double>(i)/Na);
      if (this->GetAbortExecute()) {
        break;
      }
    }

/* // @TODO fix parallel ghost cell operation

    // Check ghost points.
    // If we are processing a piece, we do not want to duplicate 
    // glyphs on the borders.  The corrct check here is:
    // ghostLevel > 0.  I am leaving this over glyphing here because
    // it make a nice example (sphereGhost.tcl) to show the 
    // point ghost levels with the glyph filter.  I am not certain 
    // of the usefullness of point ghost levels over 1, but I will have
    // to think about it.
    if (inGhostLevels && inGhostLevels[inPtId] > requestedGhostLevel) {
      continue;
    }

    if (!this->IsPointVisible(input, inPtId)) {
      continue;
      }
*/

    // Get Input point
    double *x = maskedpoints->GetPoint(i);

    // translate to Input point
    trans->Translate(x[0], x[1], x[2]);

    if (radiusdata) {
      if (!radiusMagnitude) radius = radiusdata->GetTuple1(i);
      else radius = vtkMath::Norm(radiusdata->GetTuple3(i));
    }
    if (lengthdata) {
      if (!lengthMagnitude) length = radiusdata->GetTuple1(i);
      else length = vtkMath::Norm(radiusdata->GetTuple3(i));
    }
    if (orientdata) {
      orientvector = orientdata->GetTuple3(i);
    }
  
    this->ArrowSourceObject->SetShaftRadius(radius*this->RadiusFactor);
  }
/*
  int requestedGhostLevel;
  unsigned char* inGhostLevels=0;
  vtkIdType numPts, numSourcePts, numSourceCells, inPtId, i;
  int index;
  vtkPoints *sourcePts = NULL;
  vtkSmartPointer<vtkPoints> transformedSourcePts = vtkSmartPointer<vtkPoints>::New();
  vtkPoints *newPts;
  vtkDataArray *newScalars=NULL;
  vtkDataArray *newVectors=NULL;
  vtkDataArray *newNormals=NULL;
  vtkDataArray *newTCoords = NULL;
  double x[3], v[3], vNew[3], s = 0.0, vMag = 0.0, value, tc[3];
  vtkCell *cell;
  vtkIdList *cellPts;
  int npts;
  vtkIdList *pts;
  vtkIdType ptIncr, cellIncr, cellId;
  int haveVectors, haveNormals;
  double scalex,scaley,scalez, den;
  vtkPointData* outputPD = output->GetPointData();
  vtkCellData* outputCD = output->GetCellData();
  int numberOfSources = this->GetNumberOfInputConnections(1);
  vtkPolyData *defaultSource = NULL;
  vtkIdTypeArray *pointIds=0;
  vtkPolyData *source = 0;

  vtkDebugMacro(<<"Generating glyphs");

  pts = vtkIdList::New();
  pts->Allocate(VTK_CELL_SIZE);

  pd = input->GetPointData();
  inSScalars = this->GetInputArrayToProcess(0,inputVector);
  inVectors = this->GetInputArrayToProcess(1,inputVector);
  inNormals = this->GetInputArrayToProcess(2,inputVector);
  inCScalars = this->GetInputArrayToProcess(3,inputVector);
  if (inCScalars == NULL)
    {
    inCScalars = inSScalars;
    }
  
  vtkDataArray* temp = 0;
  if (pd)
    {
    temp = pd->GetArray("vtkGhostLevels");
    }
  if ( (!temp) || (temp->GetDataType() != VTK_UNSIGNED_CHAR)
    || (temp->GetNumberOfComponents() != 1))
    {
    vtkDebugMacro("No appropriate ghost levels field available.");
    }
  else
    {
    inGhostLevels =static_cast<vtkUnsignedCharArray *>(temp)->GetPointer(0);
    }

  requestedGhostLevel =
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  
  numPts = input->GetNumberOfPoints();
  if (numPts < 1)
    {
    vtkDebugMacro(<<"No points to glyph!");
    pts->Delete();
    trans->Delete();
    return 1;
    }

  // Check input for consistency
  //
  if ( (den = this->Range[1] - this->Range[0]) == 0.0 )
    {
    den = 1.0;
    }
  if (inVectors != NULL)
    {
    haveVectors = 1;
    }
  else
    {
    haveVectors = 0;
    }

  // Allocate storage for output PolyData
  //
  outputPD->CopyVectorsOff();
  outputPD->CopyNormalsOff();
  outputPD->CopyTCoordsOff();
 
  //
  this->ArrowSource->Update();
  //

    {
    source = this->ArrowSource->GetOutput();
    sourcePts = source->GetPoints();
    numSourcePts = sourcePts->GetNumberOfPoints();
    numSourceCells = source->GetNumberOfCells();

    sourceNormals = source->GetPointData()->GetNormals();
    if ( sourceNormals )
      {
      haveNormals = 1;
      }
    else
      {
      haveNormals = 0;
      }
    
    // Prepare to copy output.
    pd = input->GetPointData();
    outputPD->CopyAllocate(pd,numPts*numSourcePts);
    if (this->FillCellData)
      {
      outputCD->CopyAllocate(pd,numPts*numSourceCells);
      }
    }

  newPts = vtkPoints::New();
  newPts->Allocate(numPts*numSourcePts);
  if ( this->GeneratePointIds )
    {
    pointIds = vtkIdTypeArray::New();
    pointIds->SetName(this->PointIdsName);
    pointIds->Allocate(numPts*numSourcePts);
    outputPD->AddArray(pointIds);
    pointIds->Delete();
    }
  if ( this->ColorMode == VTK_COLOR_BY_SCALAR && inCScalars )
    {
    newScalars = inCScalars->NewInstance();
    newScalars->SetNumberOfComponents(inCScalars->GetNumberOfComponents());
    newScalars->Allocate(inCScalars->GetNumberOfComponents()*numPts*numSourcePts);
    newScalars->SetName(inCScalars->GetName());
    }
  else if ( (this->ColorMode == VTK_COLOR_BY_SCALE) && inSScalars)
    {
    newScalars = vtkFloatArray::New();
    newScalars->Allocate(numPts*numSourcePts);
    newScalars->SetName("GlyphScale");
    if (this->ScaleMode == VTK_SCALE_BY_SCALAR)
      {
      newScalars->SetName(inSScalars->GetName());
      }
    }
  else if ( (this->ColorMode == VTK_COLOR_BY_VECTOR) && haveVectors)
    {
    newScalars = vtkFloatArray::New();
    newScalars->Allocate(numPts*numSourcePts);
    newScalars->SetName("VectorMagnitude");
    }
  if ( haveVectors )
    {
    newVectors = vtkFloatArray::New();
    newVectors->SetNumberOfComponents(3);
    newVectors->Allocate(3*numPts*numSourcePts);
    newVectors->SetName("GlyphVector");
    }
  if ( haveNormals )
    {
    newNormals = vtkFloatArray::New();
    newNormals->SetNumberOfComponents(3);
    newNormals->Allocate(3*numPts*numSourcePts);
    newNormals->SetName("Normals");
    }
    
  // Setting up for calls to PolyData::InsertNextCell()
  // JB : Fix this
    {
    output->Allocate(this->GetSource(0, inputVector[1]),
                     3*numPts*numSourceCells, numPts*numSourceCells);
    }

  transformedSourcePts->SetDataTypeToDouble();
  transformedSourcePts->Allocate(numSourcePts);

  // Traverse all Input points, transforming Source points and copying
  // point attributes.
  //
  ptIncr=0;
  cellIncr=0;
  for (inPtId=0; inPtId < numPts; inPtId++)
    {
    scalex = scaley = scalez = 1.0;
    if ( ! (inPtId % 10000) )
      {
      this->UpdateProgress(static_cast<double>(inPtId)/numPts);
      if (this->GetAbortExecute())
        {
        break;
        }
      }

    // Get the scalar and vector data
    if ( inSScalars )
      {
      s = inSScalars->GetComponent(inPtId, 0);
      if ( this->ScaleMode == VTK_SCALE_BY_SCALAR ||
           this->ScaleMode == VTK_DATA_SCALING_OFF )
        {
        scalex = scaley = scalez = s;
        }
      }
    
    if ( haveVectors )
      {
        {
        inVectors->GetTuple(inPtId, v);
        }
      vMag = vtkMath::Norm(v);
      if ( this->ScaleMode == VTK_SCALE_BY_VECTORCOMPONENTS )
        {
        scalex = v[0];
        scaley = v[1];
        scalez = v[2];
        }
      else if ( this->ScaleMode == VTK_SCALE_BY_VECTOR )
        {
        scalex = scaley = scalez = vMag;
        }
      }
    
    // Clamp data scale if enabled
    if ( this->Clamping )
      {
      scalex = (scalex < this->Range[0] ? this->Range[0] :
                (scalex > this->Range[1] ? this->Range[1] : scalex));
      scalex = (scalex - this->Range[0]) / den;
      scaley = (scaley < this->Range[0] ? this->Range[0] :
                (scaley > this->Range[1] ? this->Range[1] : scaley));
      scaley = (scaley - this->Range[0]) / den;
      scalez = (scalez < this->Range[0] ? this->Range[0] :
                (scalez > this->Range[1] ? this->Range[1] : scalez));
      scalez = (scalez - this->Range[0]) / den;
      }
    
    // Compute index into table of glyphs
    if ( this->IndexMode == VTK_INDEXING_OFF )
      {
      index = 0;
      }
    else 
      {
      if ( this->IndexMode == VTK_INDEXING_BY_SCALAR )
        {
        value = s;
        }
      else
        {
        value = vMag;
        }
      
      index = static_cast<int>((value - this->Range[0])*numberOfSources / den);
      index = (index < 0 ? 0 :
              (index >= numberOfSources ? (numberOfSources-1) : index));
      
      source = this->GetSource(index, inputVector[1]);
      if ( source != NULL )
        {
        sourcePts = source->GetPoints();
        sourceNormals = source->GetPointData()->GetNormals();
        numSourcePts = sourcePts->GetNumberOfPoints();
        numSourceCells = source->GetNumberOfCells();
        }
      }

    // Make sure we're not indexing into empty glyph
    if ( this->GetSource(index, inputVector[1]) == NULL )
      {
      continue;
      }

    // Check ghost points.
    // If we are processing a piece, we do not want to duplicate 
    // glyphs on the borders.  The corrct check here is:
    // ghostLevel > 0.  I am leaving this over glyphing here because
    // it make a nice example (sphereGhost.tcl) to show the 
    // point ghost levels with the glyph filter.  I am not certain 
    // of the usefullness of point ghost levels over 1, but I will have
    // to think about it.
    if (inGhostLevels && inGhostLevels[inPtId] > requestedGhostLevel)
      {
      continue;
      }

    if (!this->IsPointVisible(input, inPtId))
      {
      continue;
      }
    
    // Now begin copying/transforming glyph
    trans->Identity();

    // Copy all topology (transformation independent)
    for (cellId=0; cellId < numSourceCells; cellId++)
      {
      cell = this->GetSource(index, inputVector[1])->GetCell(cellId);
      cellPts = cell->GetPointIds();
      npts = cellPts->GetNumberOfIds();
      for (pts->Reset(), i=0; i < npts; i++) 
        {
        pts->InsertId(i,cellPts->GetId(i) + ptIncr);
        }
      output->InsertNextCell(cell->GetCellType(),pts);
      }
    
    // translate Source to Input point
    input->GetPoint(inPtId, x);
    trans->Translate(x[0], x[1], x[2]);
    
    if ( haveVectors )
      {
      // Copy Input vector
      for (i=0; i < numSourcePts; i++) 
        {
        newVectors->InsertTuple(i+ptIncr, v);
        }
      if (this->Orient && (vMag > 0.0))
        {
        // if there is no y or z component
        if ( v[1] == 0.0 && v[2] == 0.0 )
          {
          if (v[0] < 0) //just flip x if we need to
            {
            trans->RotateWXYZ(180.0,0,1,0);
            }
          }
        else
          {
          vNew[0] = (v[0]+vMag) / 2.0;
          vNew[1] = v[1] / 2.0;
          vNew[2] = v[2] / 2.0;
          trans->RotateWXYZ(180.0,vNew[0],vNew[1],vNew[2]);
          }
        }
      }
        
    // determine scale factor from scalars if appropriate
    // Copy scalar value
    if (inSScalars && (this->ColorMode == VTK_COLOR_BY_SCALE))
      {
      for (i=0; i < numSourcePts; i++)
        {
        newScalars->InsertTuple(i+ptIncr, &scalex); // = scaley = scalez
        }
      }
    else if (inCScalars && (this->ColorMode == VTK_COLOR_BY_SCALAR))
      {
      for (i=0; i < numSourcePts; i++)
        {
        outputPD->CopyTuple(inCScalars, newScalars, inPtId, ptIncr+i);
        }
      }
    if (haveVectors && this->ColorMode == VTK_COLOR_BY_VECTOR)
      {
      for (i=0; i < numSourcePts; i++) 
        {
        newScalars->InsertTuple(i+ptIncr, &vMag);
        }
      }
    
    // scale data if appropriate
    if ( this->Scaling )
      {
      if ( this->ScaleMode == VTK_DATA_SCALING_OFF )
        {
        scalex = scaley = scalez = this->ScaleFactor;
        }
      else
        {
        scalex *= this->ScaleFactor;
        scaley *= this->ScaleFactor;
        scalez *= this->ScaleFactor;
        }
      
      if ( scalex == 0.0 )
        {
        scalex = 1.0e-10;
        }
      if ( scaley == 0.0 )
        {
        scaley = 1.0e-10;
        }
      if ( scalez == 0.0 )
        {
        scalez = 1.0e-10;
        }
      trans->Scale(scalex,scaley,scalez);
      }

    // multiply points and normals by resulting matrix
    if (this->SourceTransform)
      {
      transformedSourcePts->Reset();
      this->SourceTransform->TransformPoints(sourcePts, transformedSourcePts);
      trans->TransformPoints(transformedSourcePts, newPts);
      }
    else
      {
      trans->TransformPoints(sourcePts,newPts);
      }
    
    if ( haveNormals )
      {
      trans->TransformNormals(sourceNormals,newNormals);
      }
    
    // Copy point data from source (if possible)
    if ( pd ) 
      {
      for (i=0; i < numSourcePts; i++)
        {
        outputPD->CopyData(pd,inPtId,ptIncr+i);
        }
      if (this->FillCellData)
        {
        for (i=0; i < numSourceCells; i++)
          {
          outputCD->CopyData(pd,inPtId,cellIncr+i);
          }
        }
      }

    // If point ids are to be generated, do it here
    if ( this->GeneratePointIds )
      {
      for (i=0; i < numSourcePts; i++)
        {
        pointIds->InsertNextValue(inPtId);
        }
      }

    ptIncr += numSourcePts;
    cellIncr += numSourceCells;
    } 
  
  // Update ourselves and release memory
  //
  output->SetPoints(newPts);
  newPts->Delete();

  if (newScalars)
    {
    int idx = outputPD->AddArray(newScalars);
    outputPD->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    newScalars->Delete();
    }

  if (newVectors)
    {
    outputPD->SetVectors(newVectors);
    newVectors->Delete();
    }

  if (newNormals)
    {
    outputPD->SetNormals(newNormals);
    newNormals->Delete();
    }

  if (newTCoords)
    {
    outputPD->SetTCoords(newTCoords);
    newTCoords->Delete();
    }
  
  output->Squeeze();
  trans->Delete();
  pts->Delete();
*/
  return 1;
}


int vtkArrowGlyphFilter::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
              outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()));
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
              outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()));
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
              outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()));
  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  return 1;
}


//----------------------------------------------------------------------------
int vtkArrowGlyphFilter::FillInputPortInformation(int port, vtkInformation *info)
{
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
    }
  return 0;
}
