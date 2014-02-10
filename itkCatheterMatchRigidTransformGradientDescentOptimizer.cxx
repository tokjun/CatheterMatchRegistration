/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef _itkCatheterMatchRigidTransformGradientDescentOptimizer_hxx
#define _itkCatheterMatchRigidTransformGradientDescentOptimizer_hxx

#include "itkCatheterMatchRigidTransformGradientDescentOptimizer.h"
#include "vnl/vnl_quaternion.h"

namespace itk
{
/**
 * Advance one Step following the gradient direction
 */
void
CatheterMatchRigidTransformGradientDescentOptimizer
::AdvanceOneStep(void)
{
  const double direction = ( m_Maximize ) ? 1.0 : -1.0;
  const ScalesType & scales = this->GetScales();

  const unsigned int nParams =  m_CostFunction->GetNumberOfParameters();

  // Make sure the scales have been set
  if ( scales.size() != nParams )
    {
    itkExceptionMacro(<< "The size of Scales is "
                      << scales.size()
                      << ", but the NumberOfParameters is "
                      << nParams
                      << ".");
    }

  DerivativeType transformedGradient(nParams);
  for ( unsigned int i = 0; i < nParams; i++ )
    {
    transformedGradient[i] = m_Gradient[i] / scales[i];
    }

  // NOTE:
  //   currentPosition[0] : Translation along the transform axis (d)
  //   currentPosition[1] : Rotation about the transform axis (positive when clockwise) (theta)
  // See itkCatheterMatchRigidTransform class
  ParametersType currentPosition = this->GetCurrentPosition();
  ParametersType newPosition(nParams);

  // Compute new translation and angle
  for ( unsigned int j = 0; j < nParams; j ++)
    {
    newPosition[j] = currentPosition[j] + direction * m_LearningRate * transformedGradient[j];
    }

  // First invoke the event, so the current position
  // still corresponds to the metric values.
  this->InvokeEvent( IterationEvent() );

  this->SetCurrentPosition(newPosition);
}
} // end namespace itk

#endif
