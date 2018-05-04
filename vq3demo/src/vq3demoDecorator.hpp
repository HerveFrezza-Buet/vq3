/*
 *   Copyright (C) 2018,  CentraleSupelec
 *
 *   Author : Hervé Frezza-Buet
 *
 *   Contributor :
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU General Public
 *   License (GPL) as published by the Free Software Foundation; either
 *   version 3 of the License, or any later version.
 *   
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *   General Public License for more details.
 *   
 *   You should have received a copy of the GNU General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 *   Contact : herve.frezza-buet@centralesupelec.fr
 *
 */


#pragma once

#include <vq3.hpp>
#include <opencv2/opencv.hpp>


namespace vq3 {
  namespace demo {
    namespace decorator {
    
      /* ########### */
      /* #         # */
      /* # Colored # */
      /* #         # */
      /* ########### */
    
      template<typename MOTHER, typename KIND> 
      struct Colored : public MOTHER {
	using decorated_type = typename MOTHER::decorated_type;
	cv::Scalar vq3_color = cv::Scalar(0, 0, 0);
	Colored(const decorated_type& val) : MOTHER(val), vq3_color(cv::Scalar(0, 0, 0)) {}
	Colored& operator=(const decorated_type& val) {this->vq3_value = val;}
      };
    
      // When we decorate a non decorated value.
      template<typename MOTHER> 
      struct Colored<MOTHER, vq3::decorator::not_decorated> {
	using decorated_type = MOTHER;
	MOTHER vq3_value;
	cv::Scalar vq3_color = cv::Scalar(0, 0, 0);
	Colored(const decorated_type& val) : vq3_value(val), vq3_color(cv::Scalar(0, 0, 0)) {}
	Colored& operator=(const decorated_type& val) {vq3_value = val;}
      };
    
      // When we decorate a decorated type with no value.
      template<typename MOTHER> 
      struct Colored<MOTHER, vq3::decorator::unvalued_decoration> : public MOTHER  {
	using decorated_type = MOTHER;
	cv::Scalar vq3_color = cv::Scalar(0, 0, 0);
	Colored() : MOTHER(), vq3_color(cv::Scalar(0, 0, 0)) {}
      };
    
      // When we decorate void.
      template<> 
      struct Colored<void, vq3::decorator::not_decorated> {
	using decorated_type = void;
	cv::Scalar vq3_color = cv::Scalar(0, 0, 0);
	Colored() : vq3_color(cv::Scalar(0, 0, 0)) {}
      };
    
      template<typename MOTHER>
      using colored = Colored<MOTHER, typename vq3::decorator::decoration<MOTHER>::value_type>;

    }
  }
}
