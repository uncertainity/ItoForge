#include <ql/quantlib.hpp>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace QuantLib;
using namespace std;

void print_1d(vector<Real> arr){
    int n = arr.size();
    for (int i = 0; i < n; i++){
        cout << arr[i] <<",";
    }
}

// int main(){
//     Calendar cal = TARGET();
//     Date today(25,February,2026);
//     Settings::instance().evaluationDate() = today;

//     Real S0    = 100.0;
//     Real v0    = 0.04;   // initial variance
//     Real mu    = 0.05;   // physical drift you want
//     Real kappa = 1.5;
//     Real theta = 0.04;
//     Real xi    = 0.3;    // "sigma" of variance process in QL
//     Real rho   = -0.7;

//     // set r = mu, q = 0, drift = r - q = mu
//     Handle<YieldTermStructure> rTS(ext::make_shared<FlatForward>(today, mu, Actual365Fixed()));
//     Handle<YieldTermStructure> qTS(ext::make_shared<FlatForward>(today, 0.0, Actual365Fixed()));
//     Handle<Quote> s0(ext::make_shared<SimpleQuote>(S0));
//     auto process = ext::make_shared<HestonProcess>(rTS,qTS,s0,v0,kappa,theta,xi,rho,HestonProcess::QuadraticExponential);

//     Time T = 1.0;
//     Size steps = 252;
//     Time dt = T/steps;
//     Size seed = 42;
//     Size factors = 2;

//     TimeGrid grid(T,steps);
//     typedef PseudoRandom::rsg_type rsg_type;
//     rsg_type rsg = PseudoRandom::make_sequence_generator(factors * grid.size()-factors, seed);
//     MultiPathGenerator<rsg_type> gen(process,grid,rsg,false);    
//     const MultiPath path = gen.next().value; 
//     vector<Real> logS(grid.size()),v(grid.size());

//     for(Size i = 0; i < grid.size(); i++){
//         logS[i] = path[0][i];
//         v[i] = path[1][i];
//     }

//     cout<<"logS path:"<<endl;
//     print_1d(logS);
//     cout<<"\n";
//     cout<<"v path:"<<endl;
//     print_1d(v);
//     cout<<"\n";
// }

int main(){
    Calendar cal = TARGET();
    Date today(25,February,2026);
    Settings::instance().evaluationDate() = today;

    Real S0 = 100.0;
    Rate r  = 0.03;   // risk-free (Q)
    Rate q  = 0.00;   // dividend yield (Q)
    DayCounter dc = Actual365Fixed();

    Handle<Quote> spot(ext::make_shared<SimpleQuote>(S0));
    Handle<YieldTermStructure> rTS(ext::make_shared<FlatForward>(today, r, dc));
    Handle<YieldTermStructure> qTS(ext::make_shared<FlatForward>(today, q, dc));

    Real v0     = 0.04;
    Real kappaQ = 1.5;
    Real thetaQ = 0.04;
    Real xi     = 0.3;
    Real rho    = -0.7;
    auto process = ext::make_shared<HestonProcess>(rTS, qTS, spot, v0, kappaQ, thetaQ, xi, rho,
                    HestonProcess::QuadraticExponential);

    auto model = ext::make_shared<HestonModel>(process);
    auto engine =  ext::make_shared<AnalyticHestonEngine>(model);
    vector <Date> maturities = {
        cal.advance(today,1,Months),
        cal.advance(today,3,Months),
        cal.advance(today,6,Months),
        cal.advance(today,12,Months)
    };
    vector <int> strikes = {80,90,100,110,120};
    // making the strike grid 
    vector <vector <Real>> C(maturities.size(),vector<Real> (strikes.size(),0));
    for (Size iT = 0; iT < maturities.size();iT++){
        auto exercise = ext::make_shared<EuropeanExercise>(maturities[iT]);
        for (Size iK = 0; iK < strikes.size(); iK ++){
            auto payoff = ext::make_shared<PlainVanillaPayoff>(Option::Call,strikes[iK]);
            VanillaOption opt(payoff,exercise);
            opt.setPricingEngine(engine);
            C[iT][iK] = opt.NPV();
        }
    }
    cout<<"Call Price grid (rows = maturity, cols = strikes)"<<endl;
    cout << setw(20) << "T\\K";
    for (Real K : strikes) cout << setw(12) << K;
    cout << "\n";

    for (Size iT = 0; iT < maturities.size(); ++iT) {
        cout << setw(12) << maturities[iT];
        for (Size iK = 0; iK < strikes.size(); ++iK) {
            cout << setw(12) << fixed << setprecision(4) << C[iT][iK];
        }
        cout << "\n";
    }

}
