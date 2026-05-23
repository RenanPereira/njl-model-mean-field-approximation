#ifndef GENERALINTEGRANDPARAMETERS_H
#define GENERALINTEGRANDPARAMETERS_H


class GeneralIntegrandParameters
{
public:
    virtual void printIntegrandVariables()
    {
        std::cout << "The method 'printIntegrandVariables' is using the default method from the class GeneralIntegrandParameters!" << "\n";
    }

    virtual void behaviourAfterFailedIntegration()
    { 
    	std::cout << "Using the default behaviour after a failed integration: abort!" << "\n";
    	abort(); 
    };
};


#endif
