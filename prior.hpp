#pragma once
#ifndef prior_hpp
#define prior_hpp

#include <map>
#include <string>

#include "rng.hpp"

class Prior {
public:

  Prior();
  virtual ~Prior();

  virtual double sample(Rng &rng) = 0;

  virtual bool valid(double v) = 0;

  virtual double pdf(double v) = 0;

  virtual double logpdf(double v) = 0;

  static Prior *load(FILE *fp);

  typedef Prior* (*prior_reader_f)(FILE *fp);

  static bool register_prior(const char *name, prior_reader_f reader);
  
private:

  static std::map<std::string, prior_reader_f> readers;
};

//
// Abstract proposal
//
class Proposal {
public:

  Proposal(Prior &prior);
  virtual ~Proposal();

  virtual Prior *get_prior();
  
  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldv,
		       double &newv,
		       double &logpriorratio) = 0;

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldv,
			      double newv) = 0;
  
  virtual double log_proposal_ratio(Rng &rng, double temperature,
				    double oldv,
				    double newv) = 0;

  static Proposal* load(FILE *fp, Prior &prior);

  typedef Proposal* (*proposal_reader_f)(FILE *fp, Prior &prior);
  
  static bool register_proposal(const char *name, proposal_reader_f reader);

protected:

  Prior &prior;

private:

  static std::map<std::string, proposal_reader_f> readers;

};

class PriorProposal {
public:

  PriorProposal(Prior *prior, Proposal *proposal);
  ~PriorProposal();
  
  double sample(Rng &rng);

  double pdf(double v);

  double logpdf(double v);

  bool propose(Rng &rng, double temperature, double oldv, double &newv, double &logpriorratio);
  
  double log_proposal(Rng &rng, double temperature,
		      double oldv,
		      double newv);
  
  double log_proposal_ratio(Rng &rng, double temperature,
			    double oldv,
			    double newv);

  Prior *get_prior()
  {
    return prior;
  }
  
  Proposal *get_proposal()
  {
    return proposal;
  }

  static PriorProposal* load(const char *filename);

private:

  Prior *prior;
  Proposal *proposal;

};

//
// Standard Priors
//

class UniformPrior : public Prior { 
public:

  UniformPrior(double vmin, double vmax);
  ~UniformPrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual double logpdf(double v);

  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;
  
  double vmin;
  double vmax;
  
};

class CosinePrior : public Prior {
public:
  CosinePrior();
  ~CosinePrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual double logpdf(double v);

  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;
  
};
  
class GaussianPrior : public Prior {
public:

  GaussianPrior(double mu, double sigma);
  ~GaussianPrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual double logpdf(double v);

  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;

  double mu;
  double sigma;
  
};

class LogNormalPrior : public Prior {
public:
  
  LogNormalPrior(double muhat, double sigmahat);
  ~LogNormalPrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual double logpdf(double v);

  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;

  double muhat;
  double sigmahat;
  
};

//
// Truncated Jeffreys prior. Basically a Jeffreys prior with limits so that
// it can be sampled from and normalized (i.e. proper).
//
class JeffreysPrior : public Prior {
public:

  JeffreysPrior(double vmin, double vmax);
  ~JeffreysPrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual double logpdf(double v);

  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;

  double vmin;
  double vmax;
  double C;
  

};

//
// Standard Proposals
//

class PriorSampleProposal : public Proposal {
public:
  
  PriorSampleProposal(Prior &prior);
  ~PriorSampleProposal();

  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldv,
		       double &newv,
		       double &logpriorratio);

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldv,
			      double newv);
  
  virtual double log_proposal_ratio(Rng &rng, double temperature,
				    double oldv,
				    double newv);

  static Proposal *read(FILE *fp, Prior &prior);

private:

  static const bool REGISTRATION;
  
};

class GaussianProposal : public Proposal {
public:

  GaussianProposal(Prior &prior, double std);
  ~GaussianProposal();

  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldv,
		       double &newv,
		       double &logpriorratio);

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldv,
			      double newv);
  
  virtual double log_proposal_ratio(Rng &rng, double temperature,
				    double oldv,
				    double newv);

  static Proposal *read(FILE *fp, Prior &prior);

private:

  static const bool REGISTRATION;
  
  double std;

};
  

#endif // prior_hpp
