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
#ifndef __itkCatheterMatchRigidTransformGradientDescentOptimizer_h
#define __itkCatheterMatchRigidTransformGradientDescentOptimizer_h

#include "itkGradientDescentOptimizer.h"

namespace itk
{
/** \class CatheterMatchRigidTransformGradientDescentOptimizer
 * \brief Implement a gradient descent optimizer
 *
 * CatheterMatchRigidTransformGradientDescentOptimizer is an extension to the
 * simple gradient descent optimizer implmented in GradientDescentOptimizer.
 * At each iteration the current position is updated according to
 *
 * p(n+1) = p(n) + learningRate * d f(p(n)) / d p(n)
 *
 * \f[
 *        p_{n+1} = p_n
 *                + \mbox{learningRate}
 *                \, \frac{\partial f(p_n) }{\partial p_n}
 * \f]
 *
 * The learning rate is a fixed scalar defined via SetLearningRate().
 * The optimizer steps through a user defined number of iterations;
 * no convergence checking is done.
 * The first four components of p are assumed to be the four components
 * of the quaternion. After each update, the quaternion is normalized to
 * have a magnitude of one. This ensures the the transform is purely rigid.
 *
 * \sa GradientDescentOptimizer
 * \ingroup Numerics Optimizers
 * \ingroup ITKOptimizers
 */
class CatheterMatchRigidTransformGradientDescentOptimizer:
  public GradientDescentOptimizer
{
public:
  /** Standard class typedefs. */
  typedef CatheterMatchRigidTransformGradientDescentOptimizer Self;
  typedef GradientDescentOptimizer                         Superclass;
  typedef SmartPointer< Self >                             Pointer;
  typedef SmartPointer< const Self >                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CatheterMatchRigidTransformGradientDescentOptimizer,
               GradientDescentOptimizer);

  /**  Parameters type.
   *  It defines a position in the optimization search space. */
  typedef Superclass::ParametersType ParametersType;

  /** Advance one step following the gradient direction. */
  virtual void AdvanceOneStep(void);

protected:
  CatheterMatchRigidTransformGradientDescentOptimizer() {}
  virtual ~CatheterMatchRigidTransformGradientDescentOptimizer() {}

private:
  //purposely not implemented
  CatheterMatchRigidTransformGradientDescentOptimizer(const Self &);
  void operator=(const Self &);
};
} // end namespace itk

#endif