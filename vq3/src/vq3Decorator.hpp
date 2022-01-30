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

#include <vq3Utils.hpp>
#include <vq3Stats.hpp>

#include <iostream>
#include <type_traits>
#include <string>
#include <atomic>

namespace vq3 {
  namespace decorator {

    /* ########## */
    /* #        # */
    /* # SFINAE # */
    /* #        # */
    /* ########## */
    
    /**
       A type tag for non decorated values.
    */
    struct not_decorated {};
    
    /**
       A type tag for a decorated type with no values (obtained from a initial 'void' type argument).
    */
    struct unvalued_decoration {};

    /**
     * A type tag for decorated values with an actual non-void value.
     */
    struct valued_decoration{};
    
    
    template<typename T> using sfinae_test = void;
    
    template<bool IS_DECORATEDD, bool HAS_VALUE> struct decoration_kind               {using value_type = valued_decoration;  };
    template<>                                   struct decoration_kind<false, false> {using value_type = not_decorated;      };
    template<>                                   struct decoration_kind<false, true > {using value_type = not_decorated;      };
    template<>                                   struct decoration_kind<true,  false> {using value_type = unvalued_decoration;};
    
    /** These templates always compile. */
    template<typename T, typename = void> struct is_decorated          : std::false_type {};
    template<typename T>                  struct is_decorated<void, T> : std::false_type {};
    template<typename T, typename = void> struct has_value             : std::false_type {};
    template<typename T>                  struct has_value<void, T>    : std::false_type {};
    
    /** These specification are tried first, but compilation may fail if T has not the expected requirement. */
    template<typename T> struct is_decorated<T, sfinae_test<typename T::decorated_type>            > : std::true_type {};
    template<typename T> struct has_value   <T, sfinae_test<decltype(std::declval<T>().vq3_value)> > : std::true_type {};

    /** This tell which is the decoration status/tag of type T (decoration<T>::value_type) */
    template<typename T> struct decoration : decoration_kind<is_decorated<T>::value, has_value<T>::value> {};

    /** This prints the kind of decoration of a type (for debugging). Use "print<T>{}();" */ 
    template<typename KIND>
    struct print {void operator()() {std::cout << "none" << std::endl;}};

    template<>
    struct print<valued_decoration> {void operator()() {std::cout << "valued decoration" << std::endl;}};

    template<>
    struct print<unvalued_decoration> {void operator()() {std::cout << "unvalued decoration" << std::endl;}};

    template<>
    struct print<not_decorated> {void operator()() {std::cout << "not decorated" << std::endl;}};
  }

  namespace spec {


    /**
     * This type is a tutorial for your own design of decorators. With
     * no template specification, KIND is supposed to be
     * vq3::decoration::valued_decoration.
     */
    template<typename MOTHER, typename KIND> 
    struct Decorator : public MOTHER {
      /** This is mandatory for a decorator.*/
      using decorated_type = typename MOTHER::decorated_type;
      /** This is the data added by the decorator. It may not be boolean. */
      bool vq3_the_decoration = false;

      /** The default constructor is required */
      Decorator() = default;

      /** The copy constructor is required */
      Decorator(const Decorator&) = default;
      
      /** This is the expected constructor for a valued decorator */
      Decorator(const decorated_type& val) : MOTHER(val), vq3_the_decoration(false) {}
      
      /** Affectation from a value has to work */
      Decorator& operator=(const decorated_type& val) {this->vq3_value = val;}
    };

    /**
       This is the specification for decorating a non-void value which is not decorated yet.
    */
    template<typename MOTHER> 
    struct Decorator<MOTHER, vq3::decorator::not_decorated> {
      /** This is mandatory for a decorator.*/
      using decorated_type = MOTHER;
      /** The value attribute must have this name.*/
      MOTHER vq3_value;
      /** This is the data added by the decorator. It may not be boolean. */
      bool vq3_the_decoration = false;

      /** The default constructor is required */
      Decorator() = default;
      /** The copy constructor is required */
      Decorator(const Decorator&) = default;
      /** This is the expected constructor for a valued decorator **/
      Decorator(const decorated_type& val) : vq3_value(val), vq3_the_decoration(false) {}
      /** Affectation from a value has to work */
      Decorator& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    /**
       This is the specification for decorating an already decorated type, for which no value is defined.
    */
    template<typename MOTHER> 
    struct Decorator<MOTHER, vq3::decorator::unvalued_decoration> : public MOTHER {
      /** This is mandatory for a decorator.*/
      using decorated_type = MOTHER;
      /** This is the data added by the decorator. It may not be boolean. */
      bool vq3_the_decoration = false;
      /** The default constructor must be implemented since no value is present, no affectation is required */
      Decorator() : MOTHER(), vq3_the_decoration(false) {}
    };
    
    /**
       This is the specification when void is decorated.
    */
    template<> 
    struct Decorator<void, vq3::decorator::not_decorated> {
      /** This is mandatory for a decorator.*/
      using decorated_type = void;
      /** This is the data added by the decorator. It may not be boolean. */
      bool vq3_the_decoration = false;
      /** The default constructor must be implemented since no value is present, no affectation is required */
      Decorator() : vq3_the_decoration(false) {}
      /** The copy constructor is required */
      Decorator(const Decorator&) = default;
    };

    /**
       This is the user friendly decorator.
    */
    template<typename MOTHER>
    using decorator = Decorator<MOTHER, typename vq3::decorator::decoration<MOTHER>::value_type>;

  }

  
    
  namespace decorator {
    
    
    /* ######## */
    /* #      # */
    /* # None # */
    /* #      # */
    /* ######## */
    
    template<typename MOTHER, typename KIND> 
    struct None : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      None() = default;
      None(const None&) = default;
      None(const decorated_type& val) : MOTHER(val) {}
      None& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER> 
    struct None<MOTHER, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      None() = default;
      None(const None&) = default;
      None(const decorated_type& val) : vq3_value(val) {}
      None& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER> 
    struct None<MOTHER, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      None() : MOTHER() {}
      None(const None&) = default;
    };
    
    // When we decorate void.
    template<> 
    struct None<void, not_decorated> {
      using decorated_type = void;
      None() = default;
      None(const None&) = default;
    };
    
    template<typename MOTHER>
    using none = None<MOTHER, typename decoration<MOTHER>::value_type>;

    
    /* ########## */
    /* #        # */
    /* # Tagged # */
    /* #        # */
    /* ########## */
    
    template<typename MOTHER, typename KIND> 
    struct Tagged : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      bool vq3_tag = false;
      Tagged() = default;
      Tagged(const Tagged&) = default;
      Tagged(const decorated_type& val) : MOTHER(val), vq3_tag(false) {}
      Tagged& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER> 
    struct Tagged<MOTHER, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      bool vq3_tag = false;
      Tagged() = default;
      Tagged(const Tagged&) = default;
      Tagged(const decorated_type& val) : vq3_value(val), vq3_tag(false) {}
      Tagged& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER> 
    struct Tagged<MOTHER, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      bool vq3_tag = false;
      Tagged() : MOTHER(), vq3_tag(false) {}
      Tagged(const Tagged&) = default;
    };
    
    // When we decorate void.
    template<> 
    struct Tagged<void, not_decorated> {
      using decorated_type = void;
      bool vq3_tag = false;
      Tagged() : vq3_tag(false) {}
      Tagged(const Tagged&) = default;
    };
    
    template<typename MOTHER>
    using tagged = Tagged<MOTHER, typename decoration<MOTHER>::value_type>;

    
    /* ########### */
    /* #         # */
    /* # Counter # */
    /* #         # */
    /* ########### */
    
    template<typename MOTHER, typename KIND> 
    struct Counter : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      std::atomic_size_t vq3_counter = 0;
      Counter() = default;
      Counter(const Counter& other) : MOTHER(other), vq3_counter(0) {}
      Counter(const decorated_type& val) : MOTHER(val), vq3_counter(0) {}
      Counter& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER> 
    struct Counter<MOTHER, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      std::atomic_size_t vq3_counter = 0;
      Counter() = default;
      // The default definition for copy constructor fails with the
      // deleted copy contructor of atomic.
      Counter(const Counter& other) : vq3_value(other.vq3_value), vq3_counter(0){}
      Counter(const decorated_type& val) : vq3_value(val), vq3_counter(0) {}
      Counter& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER> 
    struct Counter<MOTHER, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      std::atomic_size_t vq3_counter = 0;
      Counter() : MOTHER(), vq3_counter(0) {}
      Counter(const Counter& other) : MOTHER(other), vq3_counter(0) {}
    };
    
    // When we decorate void.
    template<> 
    struct Counter<void, not_decorated> {
      using decorated_type = void;
      std::atomic_size_t vq3_counter = 0;
      Counter() : vq3_counter(0) {}
      Counter(const Counter& other) : vq3_counter(0) {}
    };
    
    template<typename MOTHER>
    using counter = Counter<MOTHER, typename decoration<MOTHER>::value_type>;
    
    /* ######## */
    /* #      # */
    /* # Cost # */
    /* #      # */
    /* ######## */
    
    template<typename MOTHER, typename COST_TYPE, typename KIND> 
    struct Cost : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      COST_TYPE vq3_cost;
      Cost() = default;
      Cost(const Cost&) = default;
      Cost(const decorated_type& val) : MOTHER(val), vq3_cost() {}
      Cost& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER, typename COST_TYPE> 
    struct Cost<MOTHER, COST_TYPE, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      COST_TYPE vq3_cost;
      Cost() = default;
      Cost(const Cost&) = default;
      Cost(const decorated_type& val) : vq3_value(val), vq3_cost() {}
      Cost& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER, typename COST_TYPE> 
    struct Cost<MOTHER, COST_TYPE, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      COST_TYPE vq3_cost;
      Cost() : MOTHER(), vq3_cost() {}
      Cost(const Cost&) = default;
    };
    
    // When we decorate void.
    template<typename COST_TYPE> 
    struct Cost<void, COST_TYPE, not_decorated> {
      using decorated_type = void;
      COST_TYPE vq3_cost;
      Cost() : vq3_cost() {}
      Cost(const Cost&) = default;
    };
    
    template<typename MOTHER, typename COST_TYPE>
    using cost = Cost<MOTHER, COST_TYPE, typename decoration<MOTHER>::value_type>;
  
    
    /* ############ */
    /* #          # */
    /* # Labelled # */
    /* #          # */
    /* ############ */
    
    template<typename MOTHER, typename KIND> 
    struct Labelled : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      unsigned int vq3_label = 0;
      Labelled() = default;
      Labelled(const Labelled&) = default;
      Labelled(const decorated_type& val) : MOTHER(val), vq3_label(0) {}
      Labelled& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER> 
    struct Labelled<MOTHER, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      unsigned int vq3_label = 0;
      Labelled() = default;
      Labelled(const Labelled&) = default;
      Labelled(const decorated_type& val) : vq3_value(val), vq3_label(0) {}
      Labelled& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER> 
    struct Labelled<MOTHER, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      unsigned int vq3_label = 0;
      Labelled() : MOTHER(), vq3_label(0) {}
      Labelled(const Labelled&) = default;
    };
    
    // When we decorate void.
    template<> 
    struct Labelled<void, not_decorated> {
      using decorated_type = void;
      unsigned int vq3_label = 0;
      Labelled() : vq3_label(0) {}
      Labelled(const Labelled&) = default;
    };
    
    template<typename MOTHER>
    using labelled = Labelled<MOTHER, typename decoration<MOTHER>::value_type>;
    
    /* ############## */
    /* #            # */
    /* # Efficiency # */
    /* #            # */
    /* ############## */
    
    
    template<typename MOTHER, typename KIND> 
    struct Efficiency : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      bool vq3_efficient = true;
      Efficiency() = default;
      Efficiency(const Efficiency&) = default;
      Efficiency(const decorated_type& val) : MOTHER(val), vq3_efficient(true) {}
      Efficiency& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER> 
    struct Efficiency<MOTHER, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      bool vq3_efficient = true;
      Efficiency() = default;
      Efficiency(const Efficiency&) = default;
      Efficiency(const decorated_type& val) : vq3_value(val), vq3_efficient(true) {}
      Efficiency& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER> 
    struct Efficiency<MOTHER, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      bool vq3_efficient = true;
      Efficiency() : MOTHER(), vq3_efficient(true) {}
      Efficiency(const Efficiency&) = default;
    };
    
    // When we decorate void.
    template<> 
    struct Efficiency<void, not_decorated> {
      using decorated_type = void;
      bool vq3_efficient = true;
      Efficiency() : vq3_efficient(true) {}
      Efficiency(const Efficiency&) = default;
    };
    
    template<typename MOTHER>
    using efficiency = Efficiency<MOTHER, typename decoration<MOTHER>::value_type>;
    
    
    /* ######## */
    /* #      # */
    /* # Text # */
    /* #      # */
    /* ######## */
    
    
    template<typename MOTHER, typename KIND> 
    struct Text : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      std::string vq3_text;
      Text() = default;
      Text(const Text&) = default;
      Text(const decorated_type& val) : MOTHER(val), vq3_text() {}
      Text& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER> 
    struct Text<MOTHER, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      std::string vq3_text;
      Text() = default;
      Text(const Text&) = default;
      Text(const decorated_type& val) : vq3_value(val), vq3_text() {}
      Text& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER> 
    struct Text<MOTHER, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      std::string vq3_text;
      Text() : MOTHER(), vq3_text() {}
      Text(const Text&) = default;
    };
    
    // When we decorate void.
    template<> 
    struct Text<void, not_decorated> {
      using decorated_type = void;
      std::string vq3_text;
      Text() : vq3_text() {}
      Text(const Text&) = default;
    };
    
    template<typename MOTHER>
    using text = Text<MOTHER, typename decoration<MOTHER>::value_type>;

    
    /* ####### */
    /* #     # */
    /* # Sum # */
    /* #     # */
    /* ####### */
    
    template<typename MOTHER, typename INCREMENTABLE, typename KIND> 
    struct Sum : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      vq3::utils::accum<INCREMENTABLE, double> vq3_sum;
      Sum() = default;
      Sum(const Sum&) = default;
      Sum(const decorated_type& val) : MOTHER(val), vq3_sum() {}
      Sum& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER, typename INCREMENTABLE> 
    struct Sum<MOTHER, INCREMENTABLE, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      vq3::utils::accum<INCREMENTABLE, double> vq3_sum;
      Sum() = default;
      Sum(const Sum&) = default;
      Sum(const decorated_type& val) : vq3_value(val), vq3_sum() {}
      Sum& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER, typename INCREMENTABLE> 
    struct Sum<MOTHER, INCREMENTABLE, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      vq3::utils::accum<INCREMENTABLE, double> vq3_sum;
      Sum() : MOTHER(), vq3_sum() {}
      Sum(const Sum&) = default;
    };
    
    // When we decorate void.
    template<typename INCREMENTABLE> 
    struct Sum<void, INCREMENTABLE, not_decorated> {
      using decorated_type = void;
      vq3::utils::accum<INCREMENTABLE, double> vq3_sum;
      Sum() : vq3_sum() {}
      Sum(const Sum&) = default;
    };
    
    template<typename MOTHER, typename INCREMENTABLE>
    using sum = Sum<MOTHER, INCREMENTABLE, typename decoration<MOTHER>::value_type>;

    
    /* ########### */
    /* #         # */
    /* # GridPos # */
    /* #         # */
    /* ########### */
    
    template<typename MOTHER, typename KIND> 
    struct GridPos : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      std::pair<unsigned int, unsigned int> vq3_gridpos = {0, 0};
      GridPos() = default;
      GridPos(const GridPos&) = default;
      GridPos(const decorated_type& val) : MOTHER(val), vq3_gridpos({0, 0}) {}
      GridPos& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER> 
    struct GridPos<MOTHER, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      std::pair<unsigned int, unsigned int> vq3_gridpos = {0, 0};
      GridPos() = default;
      GridPos(const GridPos&) = default;
      GridPos(const decorated_type& val) : vq3_value(val), vq3_gridpos({0, 0}) {}
      GridPos& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER> 
    struct GridPos<MOTHER, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      std::pair<unsigned int, unsigned int> vq3_gridpos = {0, 0};
      GridPos() : MOTHER(), vq3_gridpos({0, 0}) {}
      GridPos(const GridPos&) = default;
    };
    
    // When we decorate void.
    template<> 
    struct GridPos<void, not_decorated> {
      using decorated_type = void;
      std::pair<unsigned int, unsigned int> vq3_gridpos = {0, 0};
      GridPos() : vq3_gridpos({0, 0}) {}
      GridPos(const GridPos&) = default;
    };
    
    template<typename MOTHER>
    using grid_pos = GridPos<MOTHER, typename decoration<MOTHER>::value_type>;



    
    /* ############ */
    /* #          # */
    /* # Smoother # */
    /* #          # */
    /* ############ */

    
    
    template<typename MOTHER, typename VALUE, unsigned int SAVITZKY_GOLAY_ORDER, unsigned int SAVITZKY_GOLAY_WINDOW_SIZE, unsigned int SAVITZKY_GOLAY_DEGREE, typename KIND> 
    struct Smoother : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      vq3::utils::savitzky_golay::constant_timestep::estimator<VALUE, SAVITZKY_GOLAY_ORDER, SAVITZKY_GOLAY_WINDOW_SIZE, SAVITZKY_GOLAY_DEGREE> vq3_smoother;
      Smoother() = default;
      Smoother(const Smoother&) = default;
      Smoother(const decorated_type& val) : MOTHER(val), vq3_smoother() {}
      Smoother& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER, typename VALUE, unsigned int SAVITZKY_GOLAY_ORDER, unsigned int SAVITZKY_GOLAY_WINDOW_SIZE, unsigned int SAVITZKY_GOLAY_DEGREE> 
    struct Smoother<MOTHER, VALUE, SAVITZKY_GOLAY_ORDER, SAVITZKY_GOLAY_WINDOW_SIZE, SAVITZKY_GOLAY_DEGREE, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      vq3::utils::savitzky_golay::constant_timestep::estimator<VALUE, SAVITZKY_GOLAY_ORDER, SAVITZKY_GOLAY_WINDOW_SIZE, SAVITZKY_GOLAY_DEGREE> vq3_smoother;
      Smoother() = default;
      Smoother(const Smoother&) = default;
      Smoother(const decorated_type& val) : vq3_value(val), vq3_smoother() {}
      Smoother& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER, typename VALUE, unsigned int SAVITZKY_GOLAY_ORDER, unsigned int SAVITZKY_GOLAY_WINDOW_SIZE, unsigned int SAVITZKY_GOLAY_DEGREE> 
    struct Smoother<MOTHER, VALUE, SAVITZKY_GOLAY_ORDER, SAVITZKY_GOLAY_WINDOW_SIZE, SAVITZKY_GOLAY_DEGREE, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      vq3::utils::savitzky_golay::constant_timestep::estimator<VALUE, SAVITZKY_GOLAY_ORDER, SAVITZKY_GOLAY_WINDOW_SIZE, SAVITZKY_GOLAY_DEGREE> vq3_smoother;
      Smoother() : MOTHER(), vq3_smoother() {}
      Smoother(const Smoother&) = default;
    };
    
    // When we decorate void.
    template<typename VALUE, unsigned int SAVITZKY_GOLAY_ORDER, unsigned int SAVITZKY_GOLAY_WINDOW_SIZE, unsigned int SAVITZKY_GOLAY_DEGREE> 
    struct Smoother<void, VALUE, SAVITZKY_GOLAY_ORDER, SAVITZKY_GOLAY_WINDOW_SIZE, SAVITZKY_GOLAY_DEGREE, not_decorated> {
      using decorated_type = void;
      vq3::utils::savitzky_golay::constant_timestep::estimator<VALUE, SAVITZKY_GOLAY_ORDER, SAVITZKY_GOLAY_WINDOW_SIZE, SAVITZKY_GOLAY_DEGREE> vq3_smoother;
      Smoother() : vq3_smoother() {}
      Smoother(const Smoother&) = default;
    };
    
    template<typename MOTHER, typename VALUE, unsigned int SAVITZKY_GOLAY_ORDER, unsigned int SAVITZKY_GOLAY_WINDOW_SIZE, unsigned int SAVITZKY_GOLAY_DEGREE>
    using smoother = Smoother<MOTHER, VALUE, SAVITZKY_GOLAY_ORDER, SAVITZKY_GOLAY_WINDOW_SIZE, SAVITZKY_GOLAY_DEGREE, typename decoration<MOTHER>::value_type>;

    
    /* ########## */
    /* #        # */
    /* # Custom # */
    /* #        # */
    /* ########## */
    
    template<typename MOTHER, typename CUSTOM_TYPE, typename KIND> 
    struct Custom : public MOTHER {
      using decorated_type = typename MOTHER::decorated_type;
      CUSTOM_TYPE vq3_custom;
      Custom() = default;
      Custom(const Custom&) = default;
      Custom(const decorated_type& val) : MOTHER(val), vq3_custom() {}
      Custom& operator=(const decorated_type& val) {this->vq3_value = val;}
    };
    
    // When we decorate a non decorated value.
    template<typename MOTHER, typename CUSTOM_TYPE> 
    struct Custom<MOTHER, CUSTOM_TYPE, not_decorated> {
      using decorated_type = MOTHER;
      MOTHER vq3_value;
      CUSTOM_TYPE vq3_custom;
      Custom() = default;
      Custom(const Custom&) = default;
      Custom(const decorated_type& val) : vq3_value(val), vq3_custom() {}
      Custom& operator=(const decorated_type& val) {vq3_value = val;}
    };
    
    // When we decorate a decorated type with no value.
    template<typename MOTHER, typename CUSTOM_TYPE> 
    struct Custom<MOTHER, CUSTOM_TYPE, unvalued_decoration> : public MOTHER {
      using decorated_type = MOTHER;
      CUSTOM_TYPE vq3_custom;
      Custom() : MOTHER(), vq3_custom() {}
      Custom(const Custom&) = default;
    };
    
    // When we decorate void.
    template<typename CUSTOM_TYPE> 
    struct Custom<void, CUSTOM_TYPE, not_decorated> {
      using decorated_type = void;
      CUSTOM_TYPE vq3_custom;
      Custom() : vq3_custom() {}
      Custom(const Custom&) = default;
    };
    
    template<typename MOTHER, typename CUSTOM_TYPE>
    using custom = Custom<MOTHER, CUSTOM_TYPE, typename decoration<MOTHER>::value_type>;
  }

  namespace decorator {
    namespace online {
    
      /* ##################### */
      /* #                   # */
      /* # Mean/Std (online) # */
      /* #                   # */
      /* ##################### */
    
      template<typename MOTHER, typename VALUE, typename ONLINE_PARAM, typename KIND> 
      struct MeanStd : public MOTHER {
	using decorated_type = typename MOTHER::decorated_type;
	vq3::stats::online::MeanStd<VALUE, ONLINE_PARAM> vq3_online_mean_std;
	MeanStd() = default;
	MeanStd(const MeanStd&) = default;
	MeanStd(const decorated_type& val) : MOTHER(val), vq3_online_mean_std() {}
	MeanStd& operator=(const decorated_type& val) {this->vq3_value = val;}
      };
    
      // When we decorate a non decorated value.
      template<typename MOTHER, typename VALUE, typename ONLINE_PARAM> 
      struct MeanStd<MOTHER, VALUE, ONLINE_PARAM, vq3::decorator::not_decorated> {
	using decorated_type = MOTHER;
	MOTHER vq3_value;
	vq3::stats::online::MeanStd<VALUE, ONLINE_PARAM> vq3_online_mean_std;
	MeanStd() = default;
	MeanStd(const MeanStd&) = default;
	MeanStd(const decorated_type& val) : vq3_value(val), vq3_online_mean_std() {}
	MeanStd& operator=(const decorated_type& val) {vq3_value = val;}
      };
    
      // When we decorate a decorated type with no value.
      template<typename MOTHER, typename VALUE, typename ONLINE_PARAM> 
      struct MeanStd<MOTHER, VALUE, ONLINE_PARAM, vq3::decorator::unvalued_decoration> : public MOTHER {
	using decorated_type = MOTHER;
	vq3::stats::online::MeanStd<VALUE, ONLINE_PARAM> vq3_online_mean_std;
	MeanStd() : MOTHER(), vq3_online_mean_std() {}
	MeanStd(const MeanStd&) = default;
      };
    
      // When we decorate void.
      template<typename VALUE, typename ONLINE_PARAM> 
      struct MeanStd<void, VALUE, ONLINE_PARAM, vq3::decorator::not_decorated> {
	using decorated_type = void;
	vq3::stats::online::MeanStd<VALUE, ONLINE_PARAM> vq3_online_mean_std;
	MeanStd() : vq3_online_mean_std() {}
	MeanStd(const MeanStd&) = default;
      };
    
      template<typename MOTHER, typename VALUE, typename ONLINE_PARAM>
      using mean_std = MeanStd<MOTHER, VALUE, ONLINE_PARAM, typename vq3::decorator::decoration<MOTHER>::value_type>;
    }
  }
}
