/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkArrowGlyphFilter.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkArrowGlyphFilter - A Glyph filter for arrows
// .SECTION Description
// vtkArrowGlyphFilter glyphs arrows using the independent variables
// for radius and length

#ifndef __vtkArrowGlyphFilter_h
#define __vtkArrowGlyphFilter_h

#include "vtkSmartPointer.h"
#include "vtkPolyDataAlgorithm.h"

class vtkArrowSource;
class vtkMaskPoints;

class VTK_EXPORT vtkArrowGlyphFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkArrowGlyphFilter *New();
  vtkTypeMacro(vtkArrowGlyphFilter, vtkPolyDataAlgorithm);

  // Description:
  // Turn on/off orienting of arrows along vector/normal.
  vtkSetMacro(Orient,int);
  vtkBooleanMacro(Orient,int);
  vtkGetMacro(Orient,int);

  // Description:
  // Array (vector) to use to control Orientation
  vtkSetStringMacro(OrientArray);
  vtkGetStringMacro(OrientArray);

  // Description:
  // Turn on/off scaling of arrows according to LengthArray
  vtkSetMacro(LengthScaling,int);
  vtkBooleanMacro(LengthScaling,int);
  vtkGetMacro(LengthScaling,int);

  // Description:
  // A Scaling factor to apply to the arrows in conjunction with LengthArray
  vtkSetMacro(LengthFactor,double);
  vtkGetMacro(LengthFactor,double);

  // Description:
  // Array to use to control Length
  vtkSetStringMacro(LengthArray);
  vtkGetStringMacro(LengthArray);

  // Description:
  // Turn on/off scaling of arrows according to RadiusArray
  vtkSetMacro(RadiusScaling,int);
  vtkBooleanMacro(RadiusScaling,int);
  vtkGetMacro(RadiusScaling,int);

  // Description:
  // A Scaling factor to apply to the arrows in conjunction with RadiusArray
  vtkSetMacro(RadiusFactor,double);
  vtkGetMacro(RadiusFactor,double);

  // Description:
  // Array to use to control Radius
  vtkSetStringMacro(RadiusArray);
  vtkGetStringMacro(RadiusArray);

  // Description:
  // Limit the number of points to glyph
  vtkSetMacro(MaximumNumberOfPoints, int);
  vtkGetMacro(MaximumNumberOfPoints, int);
  
  // Description:
  // Set/get whether to mask points
  void SetUseMaskPoints(int useMaskPoints);
  vtkGetMacro(UseMaskPoints, int);

  // Description:
  // Set/get flag to cause randomization of which points to mask.
  void SetRandomMode(int mode);
  int GetRandomMode();


  // Description:
  // Enable/disable the generation of point ids as part of the output. The
  // point ids are the id of the input generating point. The point ids are
  // stored in the output point field data and named "InputPointIds". Point
  // generation is useful for debugging and pick operations.
  vtkSetMacro(GeneratePointIds,int);
  vtkGetMacro(GeneratePointIds,int);
  vtkBooleanMacro(GeneratePointIds,int);

  // Description:
  // Set/Get the name of the PointIds array if generated. By default the Ids
  // are named "InputPointIds", but this can be changed with this function.
  vtkSetStringMacro(PointIdsName);
  vtkGetStringMacro(PointIdsName);

  // Description:
  // This can be overwritten by subclass to return 0 when a point is
  // blanked. Default implementation is to always return 1;
  virtual int IsPointVisible(vtkDataSet*, vtkIdType) {return 1;};

  virtual void SetArrowSourceObject(vtkArrowSource* arrowsource);

  // Description:
  // Overridden to include ArrowSourceObject's MTime.
  virtual unsigned long GetMTime();

protected:
  vtkArrowGlyphFilter();
  ~vtkArrowGlyphFilter();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int FillInputPortInformation(int, vtkInformation *);

  vtkIdType GatherTotalNumberOfPoints(vtkIdType localNumPts);
  int MaskAndExecute(vtkIdType numPts, vtkIdType maxNumPts,
                     vtkDataSet* input,
                     vtkInformation* request,
                     vtkInformationVector** inputVector,
                     vtkInformationVector* outputVector);

  int             Orient;
  char           *OrientArray;
  //
  int             LengthScaling;
  double          LengthFactor;
  char           *LengthArray;
  //
  int             RadiusScaling;
  double          RadiusFactor;
  char           *RadiusArray;
  //
  vtkMaskPoints  *MaskPoints;
  int             MaximumNumberOfPoints;
  int             UseMaskPoints;
  int             RandomMode;

  // produce input points ids for each output point
  int             GeneratePointIds; 
  char           *PointIdsName;
  
  //
//  vtkSmartPointer<vtkArrowSource> ArrowSourceObject;
  vtkArrowSource *ArrowSourceObject;

private:
  vtkArrowGlyphFilter(const vtkArrowGlyphFilter&);  // Not implemented.
  void operator=(const vtkArrowGlyphFilter&);  // Not implemented.
};

#endif
