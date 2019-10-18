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

#include <optional>
#include <chrono>

namespace vq3 {
  namespace temporal {
    class dt_averager {
    private :
      std::optional<double> dt;
      std::optional<std::chrono::time_point<std::chrono::high_resolution_clock>> time_point;
      
    public:
      
      double alpha = .1; //!< the averaging coefficient dt_ += alpha*(dt - dt_).

      /**
       * @param alpha the averaging coefficient dt_ += alpha*(dt - dt_).
       */
      dt_averager(double alpha) : dt(0), time_point(), alpha(alpha) {
	dt.reset();
	/* 
	   I initialize the attribute with dt(0), and then reset it. I
	   should have rather initialize with dt() only... I do this
	   since the default initialization triggers a "may be used
	   uninitialized" warning.
	*/
      }
      
      dt_averager()                              = default;
      dt_averager(const dt_averager&)            = default;
      dt_averager& operator=(const dt_averager&) = default;

      /**
       * Signals a time tick.
       */
      void tick() {
      	auto now      = std::chrono::high_resolution_clock::now();
      	if(time_point) {
      	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now - time_point.value()).count();
      	  if(dt)
      	    dt = *dt + alpha*(1e-6*duration - *dt);
      	  else
      	    dt = 1e-6*duration;
      	}
      	time_point = now;
      }

      void clear() {
	dt.reset();
	time_point.reset();
      }
      
      /**
       * Gives the averaged time (in seconds) between two consecutive ticks.
       */
      const std::optional<double>& operator()() const {return dt;}
	
    };
  }
}
