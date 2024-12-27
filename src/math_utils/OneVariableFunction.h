#ifndef ONEVARIABLEFUNCTION_H
#define ONEVARIABLEFUNCTION_H

class OneVariableFunction
{
private:
	double (* function)(double, void * parameters);
	void * parameters;

public:
	//OneVariableFunction(){};
	OneVariableFunction(double (* functionAux)(double, void * ), void * parametersAux){ function = functionAux; parameters = parametersAux; }
	void setFunction(double (* functionAux)(double, void * )){ function = functionAux; };
	void setParameters(void * parametersAux){ parameters = parametersAux; };
	double evaluate(double x){ return function(x, parameters); }
};


#endif