/*
    FitOCTLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FitOCTLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FitOCTLib.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.18.0

#include <stan/model/model_header.hpp>

namespace model_modHetero_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_modHetero");
    reader.add_event(22, 20, "end", "model_modHetero");
    return reader;
}

#include <meta_header.hpp>
 class model_modHetero : public prob_grad {
private:
    int N;
    vector_d x;
    vector_d y;
public:
    model_modHetero(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_modHetero(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;

        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_modHetero_namespace::model_modHetero";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 6;
            validate_non_negative_index("x", "N", N);
            context__.validate_dims("data initialization", "x", "vector_d", context__.to_vec(N));
            validate_non_negative_index("x", "N", N);
            x = vector_d(static_cast<Eigen::VectorXd::Index>(N));
            vals_r__ = context__.vals_r("x");
            pos__ = 0;
            size_t x_i_vec_lim__ = N;
            for (size_t i_vec__ = 0; i_vec__ < x_i_vec_lim__; ++i_vec__) {
                x[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(N));
            validate_non_negative_index("y", "N", N);
            y = vector_d(static_cast<Eigen::VectorXd::Index>(N));
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_i_vec_lim__ = N;
            for (size_t i_vec__ = 0; i_vec__ < y_i_vec_lim__; ++i_vec__) {
                y[i_vec__] = vals_r__[pos__++];
            }

            // validate, data variables
            current_statement_begin__ = 5;
            check_greater_or_equal(function__,"N",N,1);
            current_statement_begin__ = 6;
            current_statement_begin__ = 7;
            // initialize data variables


            // validate transformed data

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 10;
            validate_non_negative_index("theta", "2", 2);
            num_params_r__ += 2;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_modHetero() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("theta")))
            throw std::runtime_error("variable theta missing");
        vals_r__ = context__.vals_r("theta");
        pos__ = 0U;
        validate_non_negative_index("theta", "2", 2);
        context__.validate_dims("initialization", "theta", "vector_d", context__.to_vec(2));
        vector_d theta(static_cast<Eigen::VectorXd::Index>(2));
        for (int j1__ = 0U; j1__ < 2; ++j1__)
            theta(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lub_unconstrain(0,10000.0,theta);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable theta: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        typedef T__ local_scalar_t__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);

            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  theta;
            (void) theta;  // dummy to suppress unused var warning
            if (jacobian__)
                theta = in__.vector_lub_constrain(0,10000.0,2,lp__);
            else
                theta = in__.vector_lub_constrain(0,10000.0,2);


            // transformed parameters
            current_statement_begin__ = 13;
            validate_non_negative_index("sig", "N", N);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  sig(static_cast<Eigen::VectorXd::Index>(N));
            (void) sig;  // dummy to suppress unused var warning

            stan::math::initialize(sig, DUMMY_VAR__);
            stan::math::fill(sig,DUMMY_VAR__);


            current_statement_begin__ = 14;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 15;
                stan::model::assign(sig, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            (get_base1(theta,1,"theta",1) * stan::math::exp((-(get_base1(x,n,"x",1)) / get_base1(theta,2,"theta",1)))), 
                            "assigning variable sig");
            }

            // validate transformed parameters
            for (int i0__ = 0; i0__ < N; ++i0__) {
                if (stan::math::is_uninitialized(sig(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: sig" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 13;

            // model body

            current_statement_begin__ = 18;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 19;
                lp_accum__.add(normal_log<propto__>(get_base1(y,n,"y",1), 0.0, get_base1(sig,n,"sig",1)));
            }

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("theta");
        names__.push_back("sig");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(2);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;

        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "model_modHetero_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector_d theta = in__.vector_lub_constrain(0,10000.0,2);
            for (int k_0__ = 0; k_0__ < 2; ++k_0__) {
            vars__.push_back(theta[k_0__]);
            }

        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {
            current_statement_begin__ = 13;
            validate_non_negative_index("sig", "N", N);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  sig(static_cast<Eigen::VectorXd::Index>(N));
            (void) sig;  // dummy to suppress unused var warning

            stan::math::initialize(sig, DUMMY_VAR__);
            stan::math::fill(sig,DUMMY_VAR__);


            current_statement_begin__ = 14;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 15;
                stan::model::assign(sig, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            (get_base1(theta,1,"theta",1) * stan::math::exp((-(get_base1(x,n,"x",1)) / get_base1(theta,2,"theta",1)))), 
                            "assigning variable sig");
            }

            // validate transformed parameters
            current_statement_begin__ = 13;

            // write transformed parameters
            if (include_tparams__) {
            for (int k_0__ = 0; k_0__ < N; ++k_0__) {
            vars__.push_back(sig[k_0__]);
            }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities



            // validate generated quantities

            // write generated quantities
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_modHetero";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "sig" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
        }


        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= 2; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "sig" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
        }


        if (!include_gqs__) return;
    }

}; // model

}

typedef model_modHetero_namespace::model_modHetero stan_model;


#endif
