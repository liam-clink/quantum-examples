#include <iostream>
#include <iomanip>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_hermite.h>
#include <gsl/gsl_errno.h>
#include <memory>

using std::cout;
using std::endl;

struct hermite_params
{
    size_t degree;
};

double hermite(double x, void *params)
{
    auto local_params = (struct hermite_params *)params;
    return gsl_sf_hermite_func(local_params->degree, x);
}

int main()
{
    cout << "Let's do this thing!\n"
         << std::scientific;

    // Sample hermite
    for (double x = -1; x <= 1; x += .01)
    {
        cout << gsl_sf_hermite_func(1, x) << '\n';
    }

    double result;
    double error;
    size_t subdivisions = 1000;
    auto workspace = gsl_integration_workspace_alloc(subdivisions);
    hermite_params hp = {1};
    gsl_function hermite1 = {&hermite, &hp};
    int success = gsl_integration_qagi(&hermite1, 0, 0.001, subdivisions, workspace, &result, &error);

    if (success == GSL_SUCCESS)
        cout << "Result: " << result << " Absolute Error: " << error << endl;
    else
        cout << "Integration failed\n";

    gsl_integration_workspace_free(workspace);
    return 0;
}