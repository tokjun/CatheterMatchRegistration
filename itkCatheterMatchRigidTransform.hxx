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
#ifndef __itkCatheterMatchRigidTransform_hxx
#define __itkCatheterMatchRigidTransform_hxx

#include "itkCatheterMatchRigidTransform.h"

namespace itk
{
// Constructor with default arguments
template <typename TScalar>
CatheterMatchRigidTransform<TScalar>
::CatheterMatchRigidTransform() :
  Superclass(ParametersDimension)
{
  m_Rotation = VnlQuaternionType(0, 0, 0, 1); // axis * vcl_sin(t/2),
                                              // vcl_cos(t/2)
  m_BaseRotation        = VnlQuaternionType(0, 0, 0, 1);
  m_BaseTranslation     = InputVnlVectorType(0, 0, 0); 
  m_TransformAxisOrigin = InputVnlVectorType(0, 0, 0); 
  m_TransformAxisNormal = InputVnlVectorType(0, 0, 1); 

  this->m_Parameters[0] = 0.0;
  this->m_Parameters[1] = 0.0;
  this->SetParameters(this->m_Parameters);
  
}

// Constructor with default arguments
template <typename TScalar>
CatheterMatchRigidTransform<TScalar>::CatheterMatchRigidTransform(unsigned int parametersDimension) :
  Superclass(parametersDimension)
{
  m_Rotation = VnlQuaternionType(0, 0, 0, 1); // axis * vcl_sin(t/2),
                                              // vcl_cos(t/2)
  m_BaseRotation        = VnlQuaternionType(0, 0, 0, 1);
  m_BaseTranslation     = InputVnlVectorType(0, 0, 0); 
  m_TransformAxisOrigin = InputVnlVectorType(0, 0, 0); 
  m_TransformAxisNormal = InputVnlVectorType(0, 0, 1); 

  this->m_Parameters[0] = 0.0;
  this->m_Parameters[1] = 0.0;
  this->SetParameters(this->m_Parameters);
}

// Constructor with explicit arguments
template <typename TScalar>
CatheterMatchRigidTransform<TScalar>::CatheterMatchRigidTransform(const MatrixType & matrix,
                                                                const OutputVectorType & offset) :
  Superclass(matrix, offset)
{
  this->ComputeMatrixParameters();
}

// Print self
template <typename TScalar>
void
CatheterMatchRigidTransform<TScalar>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Rotation:    " << m_Rotation    << std::endl;
}

// Set rotation
template <typename TScalar>
void
CatheterMatchRigidTransform<TScalar>::SetRotation(const VnlQuaternionType & rotation)
{
  m_Rotation        = rotation;

  this->ComputeMatrix();
}

// Set the parameters in order to fit an Identity transform
template <typename TScalar>
void
CatheterMatchRigidTransform<TScalar>::SetIdentity(void)
{
  m_Rotation = VnlQuaternionType(0, 0, 0, 1);

  m_BaseTranslation = InputVnlVectorType(0, 0, 0);
  m_BaseRotation = VnlQuaternionType(0, 0, 0, 1);
  m_TransformAxisOrigin = InputVnlVectorType(0, 0, 0);
  m_TransformAxisNormal = InputVnlVectorType(0, 0, 1);

  this->Superclass::SetIdentity();
}

// Set Parameters
template <typename TScalar>
void
CatheterMatchRigidTransform<TScalar>
::SetParameters(const ParametersType & parameters)
{

  // parameter[0] : Translation along the transform axis (d)
  // parameter[1] : Rotation about the transform axis (positive when clockwise) (theta)
  
  // Let 
  //   base rotation quaternion: q_b
  //   base translation: <t_b>
  //   transform axis vector: <v_a>
  // 
  //   q = q_1 + i * q_2 + i * q_3, j * q_4
  //      q_1 = cos (theta/2)
  //      q_2 = v_ax sin (theta/2)
  //      q_3 = v_ay sin (theta/2)
  //      q_4 = v_az sin (theta/2)
  //
  // Transformation can be described as:
  //   L(<v>) = q q_b (<v> + <t_b> - <p_b>) q*_b q* + <p_b> + d <v_a>
  //          = q q_b (<v> - <p_b>) q*_b q* + <p_b> + ( q q_b <t_b> q*_b q* + d <v_a> )
  //            -----------------------------------   ---------------------------------
  //                   Rotation about <p_b>                        offset
  //
  //          = q q_b (<v>) q*_b q* + q q_b (<t_b> - <p_b>) q*_b q* + <p_b> + d <v_a> )
  //            -------------------   --------------------------------------------------
  //               Rotation                              offset


  // Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &parameters != &(this->m_Parameters) )
    {
    this->m_Parameters = parameters;
    }
  TScalar d = parameters[0];
  TScalar theta = parameters[1];

  VnlQuaternionType q(m_TransformAxisNormal, theta);

  m_Rotation = q * m_BaseRotation;

  OutputVnlVectorType vnlTranslation 
    = m_Rotation.rotate(m_BaseTranslation - m_TransformAxisOrigin)
    + m_TransformAxisOrigin + d * m_TransformAxisNormal;

  this->ComputeMatrix();

  OutputVectorType translation;
  translation.SetVnlVector(vnlTranslation);
  this->SetVarTranslation(translation);
  this->ComputeOffset();
  
  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();
}

// Set Parameters
template <typename TScalar>
const
typename CatheterMatchRigidTransform<TScalar>::ParametersType
& CatheterMatchRigidTransform<TScalar>
::GetParameters() const
  {
    
  // TODO
  VnlQuaternionType quaternion  = this->GetRotation();
  OutputVectorType  translation = this->GetTranslation();
  
  //// Transfer the quaternion part
  //unsigned int par = 0;
  //for( unsigned int j = 0; j < 4; j++ )
  //  {
  //  this->m_Parameters[par] = quaternion[j];
  //  ++par;
  //  }
  //// Transfer the constant part
  //for( unsigned int i = 0; i < SpaceDimension; i++ )
  //  {
  //  this->m_Parameters[par] = translation[i];
  //  ++par;
  //  }

  VnlQuaternionType q;
  q = quaternion * m_BaseRotation.inverse();
  this->m_Parameters[1] = q.angle();

  OutputVnlVectorType t = translation.GetVnlVector();
  OutputVnlVectorType s = t - quaternion.rotate(m_BaseTranslation - m_TransformAxisOrigin) - m_TransformAxisOrigin;
  this->m_Parameters[0] = s.magnitude();

  return this->m_Parameters;
  }

template <typename TScalar>
void
CatheterMatchRigidTransform<TScalar>
::ComputeJacobianWithRespectToParameters(const InputPointType & p, JacobianType & jacobian) const
{
  // compute derivatives with respect to rotation
  jacobian.SetSize( 3, this->GetNumberOfLocalParameters() );
  jacobian.Fill(0.0);

  // compute derivertive of rotation with respect to rotation angnle
  InputVnlVectorType p_local = 
    m_Rotation.rotate(p.GetVnlVector() + m_BaseTranslation - m_TransformAxisOrigin);
  // NOTE: (d * m_TransformAxisNormal) is not necessary because it does not change the distance between the transformed p_local and the axis
  
  // projection of p_local onto the transform axis
  InputVnlVectorType proj = m_TransformAxisNormal * dot_product(p_local, m_TransformAxisNormal);

  // Distance between p and the transform axis
  InputVnlVectorType perp = p_local - proj;
  double d = perp.magnitude();
  
  // vector of derivertive
  InputVnlVectorType dp = vnl_cross_3d(m_TransformAxisNormal, p_local);
  dp.normalize();
  
  for( unsigned int dim = 0; dim < SpaceDimension; dim++ )
    {
    jacobian[dim][0] = m_TransformAxisNormal[dim];  // translation part
    jacobian[dim][1] = dp[dim] * d;                 // rotation part
    }
}


template <typename TScalar>
void
CatheterMatchRigidTransform<TScalar>
::SetTransformAxis(InputVnlVectorType & origin, InputVnlVectorType & normal)
{
  m_TransformAxisOrigin = origin;
  m_TransformAxisNormal = normal.normalize();
}


template <typename TScalar>
void
CatheterMatchRigidTransform<TScalar>::
SetBaseTransform(InputVnlVectorType & translation, VnlQuaternionType & rotation)
{
  m_BaseRotation = rotation;
  m_BaseTranslation = translation;
}


template <typename TScalar>
const typename CatheterMatchRigidTransform<TScalar>::InverseMatrixType
& CatheterMatchRigidTransform<TScalar>::GetInverseMatrix() const
  {
  // If the transform has been modified we recompute the inverse
  if( this->InverseMatrixIsOld() )
    {
    InverseMatrixType newMatrix;
    VnlQuaternionType conjugateRotation = m_Rotation.conjugate();
    VnlQuaternionType inverseRotation = conjugateRotation.inverse();
    newMatrix = inverseRotation.rotation_matrix_transpose();
    this->SetVarInverseMatrix(newMatrix);
    }
  return this->GetVarInverseMatrix();
  }

template <typename TScalar>
void
CatheterMatchRigidTransform<TScalar>::ComputeMatrix()
{
  VnlQuaternionType conjugateRotation = m_Rotation.conjugate();
  // this is done to compensate for the transposed representation
  // between VNL and ITK
  MatrixType newMatrix;

  newMatrix = conjugateRotation.rotation_matrix_transpose();
  this->SetVarMatrix(newMatrix);
}

template <typename TScalar>
void
CatheterMatchRigidTransform<TScalar>::ComputeMatrixParameters()
{
  VnlQuaternionType quat( this->GetMatrix().GetVnlMatrix() );

  m_Rotation = quat.conjugate();
}

} // namespace

#endif
