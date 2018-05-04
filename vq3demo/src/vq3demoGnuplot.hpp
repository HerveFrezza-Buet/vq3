
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <sstream>

#include <vq3.hpp>

namespace vq3 {
  namespace demo {
    namespace gnuplot {

      class chart {
      protected:

	std::string title;
	std::string file_prefix;
	std::string term;

	std::string xlabel = "";
	std::string ylabel = "";

	bool xrange = false;
	bool yrange = false;
	double xmin, xmax, ymin, ymax;

	
	std::ofstream file;
	std::string filename;

      public:
	
	chart(const std::string& title, const std::string& file_prefix)
	  : title(title),
	    file_prefix(file_prefix),
	    term("qt") {}
	chart(const chart&) = default;
	chart(chart&&)      = default;

	void set_term_pdf() {term = "pdf";}
	void set_term_qt()  {term = "qt";}
	
	void set_xrange(double x_min, double x_max) {
	  xrange = true;
	  xmin = x_min;
	  xmax = x_max;
	}
	
	void set_yrange(double y_min, double y_max) {
	  yrange = true;
	  ymin = y_min;
	  ymax = y_max;
	}

	void set_xlabel(const std::string& label) {xlabel = label;}
	void set_ylabel(const std::string& label) {ylabel = label;}

	void open() {
	  file.exceptions(std::ios::failbit | std::ios::badbit);
	  filename = file_prefix + ".gnuplot";
	  file.open(filename);

	  file << "set term " << term << std::endl;
	  if(term == "pdf")
	    file << "set output '" << file_prefix << ".pdf'" << std::endl;
	  file << std::endl
	       << "set title \"" << title << "\"" << std::endl;
	  if(xlabel != "") file << "set xlabel \"" << xlabel << "\"" << std::endl; 
	  if(ylabel != "") file << "set ylabel \"" << ylabel << "\"" << std::endl;
	  file << std::endl;
	  if(xrange) file << "set xrange [" << xmin << ':' << xmax << ']' << std::endl;
	  if(yrange) file << "set yrange [" << ymin << ':' << ymax << ']' << std::endl;
	}

	void close() {
	  std::cout << "Gnuplot file \"" << filename << "\" generated." << std::endl;
	  file.close();
	}
      };
      
      class curve : public chart {
      private:


	std::vector< std::pair<double, double> > points;
	
      public:

	using chart::chart;
	
	auto output_iterator() {return std::back_inserter(points);}

	void make() {
	  open();

	  file << "plot '-' using 1:2 with lines notitle" << std::endl;
	  for(auto& xy : points)
	    file << xy.first << ' ' << xy.second << std::endl;

	  close();
	}
      };
      
      class histogram : public chart, public vq3::stats::histogram {
      public:

	using chart::chart;
	

	
	histogram(const std::string& title, const std::string& file_prefix)
	  : chart(title, file_prefix), vq3::stats::histogram() {}
	histogram(const histogram&) = default;
	histogram(histogram&&)      = default;
	
	
	void make() {

	  vq3::stats::histogram::make();
	  open();

	  file << "set boxwidth " << .95*(bin_max-bin_min)/bin_nb << " absolute" << std::endl
	       << "set grid" << std::endl
	       << "set style fill solid noborder" << std::endl;
	  if(!yrange) file << "set yrange[0:*]" << std::endl;
	  file << "plot '-' using 1:2 with boxes lc rgb \"#5555FF\" notitle, \\" << std::endl
	       << "     '-' using 1:2 with boxes lc rgb \"#5555AA\" title \"SCI " << (unsigned int)(sci_conf*100+.5) << "%\"" << std::endl;
	  unsigned int i=0;
	  double coef = (bin_max - bin_min)/bin_nb;
	  for(auto nb : h)
	    file << bin_min + (i++ + .5)*coef << ' ' << nb << std::endl;
	  file << 'e' << std::endl;
	  i = 0;
	  for(auto nb : h) {
	    auto val = bin_min + (i++ + .5)*coef;
	    if(sci.first <= val && val < sci.second) 
	      file << val << ' ' << nb << std::endl;
	  }

	  close();
	}
      };
      
    }
  }
}
