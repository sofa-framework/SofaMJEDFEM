/******************************************************************************
*                             SOFA MJED FEM plugin                            *
*                                (c) 2006 Inria                               *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: see Authors.md                                                     *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define SOFA_COMPONENT_FORCEFIELD_MJEDTETRAHEDRALFORCEFIELD_CPP

#include <SofaMJEDFEM/MJEDTetrahedralForceField.inl>
#include <SofaMJEDFEM/config.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>


namespace sofa::component::forcefield
{

using namespace sofa::defaulttype;

// Register in the Factory
int MJEDTetrahedralForceFieldClass = core::RegisterObject("Tetrahedral FEM model of visco-hyperelastic material using MJED")
        .add< MJEDTetrahedralForceField<Vec3Types> >()
        ;

template class SOFA_MJED_FEM_API MJEDTetrahedralForceField<Vec3Types>;


} // namespace sofa::component::forcefield
