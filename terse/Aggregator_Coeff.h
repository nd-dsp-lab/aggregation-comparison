#ifndef AGGREGATOR_COEFF_H
#define AGGREGATOR_COEFF_H


#include <cinttypes>
#include <vector>
#include <cassert>
#include <utility>
#include <chrono>
#include <cmath>

using namespace std::chrono;
using std::vector;

#include "DiscreteLaplacian.h"

static const unsigned int log2_3 = 2;

#ifdef NTL_CLIENT

#include "Polynomial.h"
#include "scheme.h"

class Aggregator_Coeff{

private:
  Parameters * ctext_parms = nullptr;
  Parameters * plain_parms = nullptr;
  int num, den;
  unsigned int num_users;
  uint64_t * delta_mod_q = NULL;
  uint64_t * t_mod_q = NULL;
  Scheme sc;
  DiscreteLaplacian dl;
  Transition * t_to_q;
  Transition * q_to_t;
  long double _beta;

  //Initialized buffers only, not singletons
  void from_buffers(const uint64_t * ctext_parms_buf, const uint64_t * plain_parms_buf, const uint64_t * delta_mod_q_buf, 
                    const uint64_t * t_mod_q_buf, 
                    const uint64_t * t_q_int, const FP_TYPE * t_q_float, const uint64_t * q_t_int, const FP_TYPE * q_t_float){
    if(ctext_parms != nullptr){
      delete ctext_parms;
      ctext_parms = nullptr;
    }
    ctext_parms = new Parameters(ctext_parms_buf);

    if(plain_parms != nullptr){
      delete plain_parms;
      plain_parms = nullptr;
    }
    plain_parms = new Parameters(plain_parms_buf);

    if(delta_mod_q != NULL){
      free(delta_mod_q);
      delta_mod_q = NULL;
    }
    delta_mod_q = (uint64_t *) malloc(ctext_parms->moduli_count()*sizeof(uint64_t));
    memcpy(delta_mod_q, delta_mod_q_buf, ctext_parms->moduli_count()*sizeof(uint64_t));

    if(t_mod_q != NULL){
      free(t_mod_q);
      t_mod_q = NULL;
    }
    t_mod_q = (uint64_t *) malloc(ctext_parms->moduli_count()*sizeof(uint64_t));
    memcpy(t_mod_q, t_mod_q_buf, ctext_parms->moduli_count()*sizeof(uint64_t));

    if(t_to_q != nullptr){
      delete t_to_q;
      t_to_q = nullptr;
    }
    t_to_q = new Transition(t_q_int, t_q_float);

    if(q_to_t != nullptr){
      delete q_to_t;
      q_to_t = nullptr;
    }
    q_to_t = new Transition(q_t_int, q_t_float);
    
    return;
  }  


  void to_buffers(uint64_t ** ctext_parms_buf, uint64_t ** plain_parms_buf, uint64_t ** delta_mod_q_buf, uint64_t ** t_mod_q_buf, 
                  uint64_t ** t_q_int, FP_TYPE ** t_q_float, uint64_t ** q_t_int, FP_TYPE ** q_t_float, int & num_in, int & den_in,
                  unsigned int & n_users_in, long double & beta, Scheme & sc_in){
    //TODO also delete preexisting buffer
    if(*ctext_parms_buf != nullptr){
      free(*ctext_parms_buf);
      *ctext_parms_buf = nullptr;
    }
    if(*plain_parms_buf != nullptr){
      free(*plain_parms_buf);
      *plain_parms_buf = nullptr;
    }
    *ctext_parms_buf = ctext_parms->to_buffer();
    *plain_parms_buf = plain_parms->to_buffer();

    if(*delta_mod_q_buf != NULL){
      free(*delta_mod_q_buf);
      *delta_mod_q_buf = NULL;
    }
    *delta_mod_q_buf = (uint64_t *) malloc(ctext_parms->moduli_count()*sizeof(uint64_t));
    memcpy(*delta_mod_q_buf, delta_mod_q, ctext_parms->moduli_count()*sizeof(uint64_t));

    if(*t_mod_q_buf != NULL){
      free(*t_mod_q_buf);
      *t_mod_q_buf = NULL;
    }
    *t_mod_q_buf = (uint64_t *) malloc(ctext_parms->moduli_count()*sizeof(uint64_t));
    memcpy(*t_mod_q_buf, t_mod_q, ctext_parms->moduli_count()*sizeof(uint64_t));

    q_to_t->to_buffers(q_t_int, q_t_float);
    t_to_q->to_buffers(t_q_int, t_q_float);

    num_in = num;
    den_in = den;
    n_users_in = num_users;
    beta = _beta;
    sc_in = sc;

    return;
  }


  Aggregator_Coeff(const uint64_t * ctext_parms_buf, const uint64_t * plain_parms_buf, const uint64_t * delta_mod_q_buf, 
                    const uint64_t * t_mod_q_buf, 
                    const uint64_t * t_q_int, const FP_TYPE * t_q_float, const uint64_t * q_t_int, const FP_TYPE * q_t_float,
                    const int num_in, const int den_in,
                    const int n_users_in, const long double beta, const Scheme sc_in) {
    dl = DiscreteLaplacian(beta); 
    _beta = beta;
    sc = sc_in;
    num_users = n_users_in;
    num = num_in;
    den = den_in;
    this->from_buffers(ctext_parms_buf, plain_parms_buf, delta_mod_q_buf, 
                    t_mod_q_buf, 
                    t_q_int, t_q_float, q_t_int, q_t_float);
  }



  //Lower than necessary - just in case
  static const unsigned int PLAIN_MOD_SIZE_MAX = 50;

  std::pair<Parameters *, Parameters *> parms_ptrs() const {
    std::pair<Parameters *, Parameters *> ret;
    ret.first = ctext_parms;
    ret.second = plain_parms;
    return ret;
  }

  std::pair<Transition *, Transition *> trans_ptrs() const {
    std::pair<Transition *, Transition *> ret;
    ret.first = t_to_q;
    ret.second = q_to_t;
    return ret;
  }


  DiscreteLaplacian * dist() {
    return &dl;
  }

  unsigned int user_count() const {
    return num_users;
  }

  ~Aggregator_Coeff(){
    free(delta_mod_q);
    delta_mod_q = NULL;
    free(t_mod_q);
    t_mod_q = NULL;
    delete t_to_q;
    t_to_q = nullptr;
    delete q_to_t;
    q_to_t = nullptr;
    delete ctext_parms;
    ctext_parms = nullptr;
    delete plain_parms;
    plain_parms = nullptr;  
  }

  Aggregator_Coeff operator=(const Aggregator_Coeff & a) = delete;
  Aggregator_Coeff(const Aggregator_Coeff & a) = delete;
  Aggregator_Coeff() = delete;

  //Declare things with preexisting parameters - for small tests
  Aggregator_Coeff(const Parameters & parms_in, const uint64_t t_in, const float scale_in, 
                 const unsigned int n_users, const Scheme sc_in, 
                 const long double beta): dl(beta){
#if !USE_NTL
    assert("This operation not supported without NTL" && 0);
#else
    _beta = beta;
    sc = sc_in;
    float_to_frac(scale_in, num, den);
    num_users = n_users;
    ctext_parms = new Parameters(parms_in); //Still don't know why the compiler allows this
    assert(ctext_parms->no_repeating());
    ZZ q = ctext_parms->modulus();
    unsigned int N = ctext_parms->poly_mod_degree();
    plain_parms = new Parameters(N, {t_in}, false); //No need to multiply
    ZZ t = plain_parms->modulus();

    size_t k = ctext_parms->moduli_count();

    delta_mod_q = (uint64_t *) malloc(sizeof(uint64_t) * k);
    t_mod_q = (uint64_t *) malloc(sizeof(uint64_t) * k);
    ZZ tmp;
    ZZ delta = q/t;
    ZZ tmp_mod;
    //Fill delta mod q for later scaling
    for(size_t i = 0; i < k; i++){
      uint64_t qi = ctext_parms->moduli(i);
      uint64_t tmp_delta, tmp_t;
      //TODO fix - don't write directly to array
      tmp_mod = NTL::ZZFromBytes((const unsigned char *) &qi, sizeof(uint64_t));
      tmp = delta % tmp_mod;
      BytesFromZZ((unsigned char *) &tmp_delta, tmp, sizeof(uint64_t));
      delta_mod_q[i] = tmp_delta;
      tmp = t % tmp_mod;
      BytesFromZZ((unsigned char *) &tmp_t, tmp, sizeof(uint64_t));
      t_mod_q[i] = tmp_t;
    }

    //Transition data
    t_to_q = new Transition(*plain_parms, *ctext_parms);
    q_to_t = new Transition(*ctext_parms, *plain_parms);
#endif    
  }

  void secret_keys(Polynomial & agg_key, vector<Polynomial> & secret_keys){
    agg_key.zero();
    secret_keys.clear();
    secret_keys.reserve(num_users);

    Polynomial result_template(agg_key.parameters());
    secret_keys.resize(num_users, result_template);
    for(unsigned int i = 0; i < num_users; i++){
#ifdef DEBUG
      assert(secret_keys[i].buffer() != NULL);
#endif      
      secret_keys[i].error(this->dl);
      agg_key -= secret_keys[i];
    }
    return;
  }


  Polynomial public_key(const uint64_t ts){
    Polynomial pk(this->ctext_parms);
#ifdef DEBUG
    assert(pk.parameters() == this->ctext_parms);
#endif    
    dl.refresh(ts);
    pk.uniform(dl);
    return pk;
  }

  vector<Polynomial> public_keys(const uint64_t initial_ts, const unsigned int num_aggregations){
    Polynomial pk(this->ctext_parms);
#ifdef DEBUG
    assert(pk.parameters() == this->ctext_parms);
#endif   
    vector<Polynomial> ret;
    unsigned int total_polynomials = num_aggregations / this->ctext_parms->poly_mod_degree();
    if(num_aggregations % this->ctext_parms->poly_mod_degree()){
      total_polynomials++;
    }
    ret.reserve(total_polynomials);
    for(unsigned int i = 0; i < total_polynomials; i++){
      pk = this->public_key(initial_ts + i);
      ret.push_back(pk);
    }
    return ret;
  }

  //pk_sk_vec is a precomputed list of pk[i]*sk
  vector<uint64_t> enc(const uint64_t user_input, const unsigned int agg_idx, 
    const vector<Polynomial> & pk_sk_vec, bool do_noise, 
    double & noise_time, double & enc_time){
    steady_clock::time_point start, end;
    
    uint64_t noisy_input;

    if(do_noise){
    
      start = steady_clock::now();

      int e = dl.dl(num, den);
      uint64_t qi = this->plain_parms->moduli(0);
      uint64_t noise = e >= 0? e : qi-(-e);
      noisy_input = user_input + noise;
      noisy_input %= qi;
   
      end = steady_clock::now();
      noise_time = duration_cast<chrono::nanoseconds>(end-start).count();
    }
    else{
      noisy_input = user_input;
      noise_time = 0.0;
    }
    //These parts are not input-dependent
    vector<uint64_t> ret(this->ctext_parms->moduli_count());
    unsigned int pmd = this->ctext_parms->poly_mod_degree();
    unsigned int key_idx = agg_idx / pmd;
    unsigned int coeff_idx = agg_idx - key_idx;

    //Raise input to base q and add
    start = steady_clock::now();
    //Add noise
    int e = this->dl.uniform(3);
    
    for(size_t i = 0; i < ret.size(); i++){
      uint64_t qi = this->ctext_parms->moduli(i);
      //Move to base q, deferring reduction until after multiplication by delta
      ret[i] = noisy_input;
      
      if(this->sc == MS){
        ret[i] *= delta_mod_q[i];
      }

      ret[i] %= qi;

      uint64_t e_64 ;
      //Scale either message or noise
      if(this->sc == MS){
        e_64 = (e != 3) ? e : qi-1;;
      }
      else{
        if(e == 3){
          e_64 = qi-1;
        }
        else{
          e_64 = e * delta_mod_q[i];
        }
      }

      //Add in error term
      ret[i] += e_64;
      //Add in this user's key
      ret[i] += pk_sk_vec[key_idx].at(coeff_idx, i);
      //Finally, reduce!
      ret[i] %= qi;
    }
    end = steady_clock::now();
    enc_time = duration_cast<chrono::nanoseconds>(end-start).count();

    return ret;
  }

  uint64_t aggregate_and_decrypt(const vector<Polynomial> & agg_keys, vector<vector<uint64_t> > & ctexts, const unsigned int agg_idx,
      double & dec_time) const {
    vector<uint128_t> tmp(ctexts[0].size(), 0);
    //First, aggregate
    for(const vector<uint64_t> & v : ctexts){
      for(size_t i = 0; i < v.size(); i++){
        tmp[i] += v[i];
      }
    }
    for(size_t i = 0; i < tmp.size(); i++){
      tmp[i] %= this->ctext_parms->moduli(i);
    }
    return (this->sc == NS) ? dec_ns(agg_keys, tmp, agg_idx, dec_time) : dec_ms(agg_keys, tmp, agg_idx, dec_time);
  }

  //TODO add A*sk
  //Hardcoded for 1 plaintext modulus
  uint64_t dec_ns(const vector<Polynomial> & agg_keys, const vector<uint128_t> & aggregated_ctexts, const unsigned int agg_idx, 
    double & dec_time) const {
    steady_clock::time_point start, end;
    unsigned int key_idx = agg_idx / this->ctext_parms->poly_mod_degree();
    unsigned int coeff_idx = agg_idx - key_idx;
    uint64_t p = this->plain_parms->moduli(0);
    
    start = steady_clock::now();    

    size_t k = this->ctext_parms->moduli_count();
    FP_TYPE z;
    vector<uint64_t> y(this->ctext_parms->poly_mod_degree());

    for(size_t mod_idx = 0; mod_idx < k; mod_idx++){
      uint64_t q_tilde = this->t_to_q->q_tilde(mod_idx);
      uint64_t modulus = this->ctext_parms->moduli(mod_idx);

      //Do multiplication
      uint128_t tmp_128 = aggregated_ctexts[mod_idx] + agg_keys[key_idx].at(coeff_idx, mod_idx);
      tmp_128 *= q_tilde;
      tmp_128 %= modulus;
      y[mod_idx] = tmp_128;
      FP_TYPE z_val = y[mod_idx] / modulus;
      if(!mod_idx){
        z = z_val;
      } else{
        z += z_val;
      }
    }

    int v = std::round(z);
    uint64_t v_mod_q = v < 0 ? p : v;
    uint64_t neg_v_mod_q = v_mod_q ? p - v_mod_q : 0;
    uint64_t ret = 0;

    for(size_t i = 0; i < y.size(); i++){
      uint64_t q_star = this->t_to_q->q_star(i, 0);
      ret += y[i]*q_star;
    }
    ret += neg_v_mod_q;
    ret %= p;
    
    end = steady_clock::now();
    dec_time = duration_cast<chrono::nanoseconds>(end-start).count();
    return ret;
  }


  uint64_t dec_ms(const vector<Polynomial> & agg_keys, const vector<uint128_t> & aggregated_ctexts, const unsigned int agg_idx, 
      double & dec_time) const {
    steady_clock::time_point start, end;
    unsigned int key_idx = agg_idx / this->ctext_parms->poly_mod_degree();
    unsigned int coeff_idx = agg_idx - key_idx;
    unsigned int k = this->ctext_parms->poly_mod_degree();
    uint64_t p = this->plain_parms->moduli(0);

    int v = 0;
    FP_TYPE tmp = 0.0;
    uint128_t w = 0;

    start = steady_clock::now();

    for(size_t i = 0; i < k; i++){
      uint128_t input_plus_key = aggregated_ctexts[i] + agg_keys[key_idx].at(coeff_idx, i);
      tmp += input_plus_key * this->t_to_q->theta(i);
      uint128_t tmp_128 = input_plus_key * this->t_to_q->w(i, 0);
      tmp_128 %= p;
      w += tmp_128;
    }
    v = std::round(tmp);
    if(v < 0){
      v = p - (-v);
    }
    w += v;
    w %= p;
    end = steady_clock::now();
    dec_time = duration_cast<chrono::nanoseconds>(end-start).count();
    return w;
  }

};

vector<uint64_t> plain_integer_keys(const vector<Polynomial> & x){
  vector<uint64_t> ret;
  if(x.empty()){
    return ret;
  }
  ret.reserve(x.size()*x.front().poly_mod_degree()*x.front().mod_count());
  for(const Polynomial & p : x){
    size_t deg = p.poly_mod_degree();
    size_t towers = p.mod_count();
    for(size_t i = 0; i < deg; i++){
      for(size_t j = 0; j < towers; j++){
        ret.push_back(p.at(i, j));
      }
    }
  }
  return ret;
}

#endif


//Function signature from hell, but the deadline is looming
vector<uint64_t> enc_noclass(const uint64_t user_input, const unsigned int agg_idx, 
    const vector<uint64_t> & pk_sk_vec, bool do_noise, 
    double & noise_time, double & enc_time, const vector<uint64_t> & moduli, 
    const vector<uint64_t> & delta_mod_q,
    DiscreteLaplacian & dl, const uint64_t plain_modulus,
     int num, int den, bool noise_scaled = true){
    
    //These parts are not input-dependent
    steady_clock::time_point start, end;
    vector<uint64_t> ret(moduli.size());
    uint64_t noisy_input;

    if(do_noise){
    
      start = steady_clock::now();

      int dp = dl.dl(num, den);
      uint64_t noise = dp >= 0? dp : plain_modulus-(-dp);
      noisy_input = user_input + noise;
      noisy_input %= plain_modulus;
   
      end = steady_clock::now();
      noise_time = duration_cast<std::chrono::nanoseconds>(end-start).count();
    }
    else{
      noisy_input = user_input;
      noise_time = 0.0;
    }
    

    //Raise input to base q and add
    start = steady_clock::now();
    //Get noise
    int e = dl.uniform(3);
    
    for(size_t i = 0; i < ret.size(); i++){
      uint64_t qi = moduli[i];
      //Move to base q, deferring reduction until after multiplication by delta
      ret[i] = noisy_input;
      
      if(noise_scaled){
        ret[i] *= delta_mod_q[i];
      }

      ret[i] %= qi;

      uint64_t e_64 ;
      //Scale either message or noise
      if(!noise_scaled){
        e_64 = (e != 3) ? e : qi-1;;
      }
      else{
        if(e == 3){
          e_64 = qi-1;
        }
        else{
          e_64 = e * delta_mod_q[i];
        }
      }

      //Add in error term
      ret[i] += e_64;
      //Add in this user's key
      ret[i] += pk_sk_vec[agg_idx + i];
      //Finally, reduce!
      ret[i] %= qi;
    }
    end = steady_clock::now();
    enc_time = duration_cast<std::chrono::nanoseconds>(end-start).count();

    return ret;
  }

  //Warning: this uses __gcd, hope this works for cross-compilation
  void next_coprime(uint64_t & val, const vector<uint64_t> & arr){
    while(val){
      bool found_factor = false;
      for(const uint64_t x : arr){
        if(std::gcd(val, x) != 1){
          found_factor = true;
          break;
        }
      }
      if(!found_factor){
        return;
      }
      else{
        val++;
      }
    }
  }

  unsigned int generate_params(const uint64_t plain_modulus, const unsigned int num_users, const unsigned int num_aggregations,
    vector<uint64_t> & moduli, vector<uint64_t> & keys, vector<uint64_t> & delta_mod_q){
    const static constexpr unsigned int MOD_BITS = 59;
    //const static uint64_t MODULUS_BOUNDS = uint64_t(1) << MOD_BITS;
    //Choose for SLAP_NS
    unsigned int ctext_bits = ceil(log2(num_users) 
      + log2(plain_modulus) + log2_3);
    unsigned int num_moduli = ctext_bits / MOD_BITS;
    if(ctext_bits % MOD_BITS){
      num_moduli++;
    }

    DiscreteLaplacian dl;
    moduli.reserve(num_moduli);
    //Technically a coprime, not a prime
    uint64_t pr = 1;
    pr <<= MOD_BITS;
    for(size_t i = 0; i < moduli.size(); i++){
      moduli.push_back(pr);
      pr++;
      next_coprime(pr, moduli);
    }
    delta_mod_q.resize(num_moduli);
    for(size_t i = 0; i < delta_mod_q.size(); i++){
      delta_mod_q[i] = dl.uniform_64(moduli[i]);
    }
    keys.resize(num_aggregations*num_moduli);
    for(unsigned int i = 0; i < num_aggregations; i++){
      for(unsigned int j = 0; j < num_moduli; j++){
        keys.at((i*num_moduli) + j) = dl.uniform_64(moduli[i]);
      }
    }
    return num_moduli * MOD_BITS;
  }


#endif