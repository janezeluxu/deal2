/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 - 2017 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Zelu Xu, 2018
 */

#include "../include/includeall.h"
using namespace std;
namespace incompressible
{
    using namespace dealii;

    double Re;// = 400;
    int maxNewtonIter ;//= 20;
    double tolerance ;//= 1e-8;
    int max_cycle ;//= 2;
    //std::string vtkfile;// = "kellysolution-";
    std::string meshfile;// = "./mesh/lidCavity001.msh";
    bool KellyError;//=true;
    //int component_mask [4] = {3,3,4,3};
    vector<int> component_mask;//{3,3,4,3};
    vector<int> component_geotag;
    int nBC_com;// = 4;
    vector<int> BC_u;
    vector<double> BCval_u;
    vector<int> BC_v;
    vector<double> BCval_v;
    vector<int> BC_p;
    vector<double> BCval_p;
          
    
// ---------not needed for now -----------
    unsigned int dim = 2;
    int unsigned max_degree = 1;
    int unsigned min_degree = 1;
    int meshRefineLevel = 0;
    unsigned int nstp = 1;
// ---------not needed for now --------

    // save to solution filename
    std::ostringstream velocity_sol("", ios_base::app);
    std::ostringstream pressure_sol("", ios_base::app);
    std::ostringstream vtkfilename("", ios_base::app);

    std::ostringstream outputmeshfilename("", ios_base::app);
    
    // read from reference solution in getTrueError
    std::ostringstream filenameRefu("./solution/ref-un1.txt");
    std::ostringstream filenameRefp("./solution/ref-pn1.txt");

    int degreeRef = 1;
    int n_Ref = 80;
    
    
    //std::ostringstream vtkfilename("kellysolution-");
    //std::string meshfile = "./mesh/lidCavity001.msh";
    bool calcError=false;
    //bool KellyError=true;
    
        
    void read_input(std::string inputfile)
    {
        // Create a text string, which is used to output the text file
        string myText;

        // Read from the text file
        //std::string inputfile ="./input/input-file.txt";
        ifstream MyReadFile(inputfile);
        std::cout << "read input " <<std::endl;
        // Use a while loop together with the getline() function to read the file line by line
        while (getline (MyReadFile, myText)) {
          // Output the text from the file
            //std::cout << "read LINE " <<std::endl;
            //std::cout << myText<<std::endl;
            
            //std::string s = "scott>=tiger";
            std::string delimiter = "=";
            size_t pos = myText.find(delimiter);
            std::string token = myText.substr(0, pos);
            //std::cout << "token "<<token<<std::endl;
            
            myText.erase(0, pos + delimiter.length());
            
            //std::string token1 = "Re";
            //std::cout << token<<token1<<std::endl;
            if (token.compare("Re ") == 0)
            {
                Re = std::stod(myText);
                std::cout << "Re "<<Re<<std::endl;
            }
            
            if (token.compare("maxNewtonIter ") == 0)
            {
                maxNewtonIter = std::stoi(myText);
                std::cout << "maxNewtonIter "<<maxNewtonIter<<std::endl;
            }
            
            if (token.compare("tolerance ") == 0)
            {
                tolerance = std::stod(myText);
                std::cout << "tolerance "<<tolerance<<std::endl;
            }
            
            if (token.compare("max_cycle ") == 0)
            {
                max_cycle = std::stoi(myText);
                std::cout << "max_cycle "<<max_cycle<<std::endl;
            }
            
            if (token.compare("vtkfilename ") == 0)
            {
                vtkfilename<< myText;
            }
            
            if (token.compare("velocity_solname ") == 0)
            {
                velocity_sol<< myText;
            }
            
            if (token.compare("pressure_solname ") == 0)
            {
                pressure_sol<< myText;
            }
            
            if (token.compare("outputmeshfilename ") == 0)
            {
                outputmeshfilename<< myText;
            }
            
            if (token.compare("meshfile ") == 0)
            {
                meshfile = myText;//"./mesh/lidCavity001.msh";
                //std::cout << "meshfile="<<meshfile<<myText<<std::endl;
            }
            
            if (token.compare("KellyError ") == 0)
            {
                int KellyErrorInt = std::stoi(myText);
                std::cout << "KellyError "<<KellyErrorInt<<std::endl;
                if (KellyErrorInt ==1)
                    KellyError = true;
                else
                    KellyError = false;
            }
            
            if (token.compare("nBC_com ") == 0)
            {
                nBC_com = std::stoi(myText);
                std::cout << "nBC_com "<<nBC_com<<std::endl;
            }
            
            if (token.compare("component_mask ") == 0)
            {
                size_t pos = 0;
                //std::string token;
                std::string delimiter = ";";
                //vector<int> component_mask_test;
                while ((pos = myText.find(delimiter)) != std::string::npos) {
                    token = myText.substr(0, pos);
                    //std::cout << token << std::endl;
                    myText.erase(0, pos + delimiter.length());
                    component_mask.push_back(std::stoi(token));
                }
                
                //for (int i = 0; i<nBC_com;i++)
                  //  std::cout << "component_mask_test "<<component_mask_test[i]<<std::endl;
            }
            
            if (token.compare("component_geotag ") == 0)
            {
                size_t pos = 0;
                //std::string token;
                std::string delimiter = ";";
                //vector<int> component_mask_test;
                while ((pos = myText.find(delimiter)) != std::string::npos) {
                    token = myText.substr(0, pos);
                    //std::cout << token << std::endl;
                    myText.erase(0, pos + delimiter.length());
                    component_geotag.push_back(std::stoi(token));
                }
                
                for (int i = 0; i<nBC_com;i++)
                    std::cout << "component_geotag "<<component_geotag[i]<<std::endl;
            }
            
            if (token.compare("BC_u ") == 0)
            {
                size_t pos = 0;
                //std::string token;
                std::string delimiter = ";";
                //vector<int> component_mask_test;
                while ((pos = myText.find(delimiter)) != std::string::npos) {
                    token = myText.substr(0, pos);
                    //std::cout << token << std::endl;
                    myText.erase(0, pos + delimiter.length());
                    BC_u.push_back(std::stoi(token));
                }
                
                //for (int i = 0; i<nBC_com;i++)
                  //  std::cout << "BC_u "<<BC_u[i]<<std::endl;
            }
            
            if (token.compare("BCval_u ") == 0)
            {
                size_t pos = 0;
                //std::string token;
                std::string delimiter = ";";
                //vector<int> component_mask_test;
                while ((pos = myText.find(delimiter)) != std::string::npos) {
                    token = myText.substr(0, pos);
                    //std::cout << token << std::endl;
                    myText.erase(0, pos + delimiter.length());
                    BCval_u.push_back(std::stod(token));
                }
                
                //for (int i = 0; i<nBC_com;i++)
                  //  std::cout << "BCval_u "<<BCval_u[i]<<std::endl;
            }
            
            if (token.compare("BC_v ") == 0)
            {
                size_t pos = 0;
                //std::string token;
                std::string delimiter = ";";
                //vector<int> component_mask_test;
                while ((pos = myText.find(delimiter)) != std::string::npos) {
                    token = myText.substr(0, pos);
                    //std::cout << token << std::endl;
                    myText.erase(0, pos + delimiter.length());
                    BC_v.push_back(std::stoi(token));
                }
                
                //for (int i = 0; i<nBC_com;i++)
                  //  std::cout << "BC_v "<<BC_v[i]<<std::endl;
            }
            
            if (token.compare("BCval_v ") == 0)
            {
                size_t pos = 0;
                //std::string token;
                std::string delimiter = ";";
                //vector<int> component_mask_test;
                while ((pos = myText.find(delimiter)) != std::string::npos) {
                    token = myText.substr(0, pos);
                    //std::cout << token << std::endl;
                    myText.erase(0, pos + delimiter.length());
                    BCval_v.push_back(std::stod(token));
                }
                
                //for (int i = 0; i<nBC_com;i++)
                  //  std::cout << "BCval_v "<<BCval_v[i]<<std::endl;
            }
            
            if (token.compare("BC_p ") == 0)
            {
                size_t pos = 0;
                //std::string token;
                std::string delimiter = ";";
                //vector<int> component_mask_test;
                while ((pos = myText.find(delimiter)) != std::string::npos) {
                    token = myText.substr(0, pos);
                    //std::cout << token << std::endl;
                    myText.erase(0, pos + delimiter.length());
                    BC_p.push_back(std::stoi(token));
                }
                
                //for (int i = 0; i<1;i++)
                  //  std::cout << "BC_p "<<BC_p[i]<<std::endl;
            }
            
            if (token.compare("BCval_p ") == 0)
            {
                size_t pos = 0;
                //std::string token;
                std::string delimiter = ";";
                //vector<int> component_mask_test;
                while ((pos = myText.find(delimiter)) != std::string::npos) {
                    token = myText.substr(0, pos);
                    //std::cout << token << std::endl;
                    myText.erase(0, pos + delimiter.length());
                    BCval_p.push_back(std::stod(token));
                }
                
                //for (int i = 0; i<1;i++)
                  //  std::cout << "BCval_p "<<BCval_p[i]<<std::endl;
            }
            
        }

        // Close the file
        MyReadFile.close();
    }
        
    double Dirichlet_u(const dealii::Point<2> & p,int BC_ID)
    {
        //std::cout << "BC_u "<<BC_u[0]<<std::endl;
        for (int i = 0; i<BC_u.size();i++)
        {
            //if (BC_ID == 1)
            //{   //return -p[1]*(p[1]-0.2); //poiseuille
            //    return 8*1.5*(2*(p[1]-0.5)-1)*(-(p[1]-0.5)); // backstep
            //}
            //else
            if (BC_ID == BC_u[i])
                return BCval_u[i];
        }
          //std::cout << "BC_u "<<BC_u[i]<<std::endl;
        //if (BC_ID == 0)
          //  return 1.0;
        return 0;
    }
    
    double Dirichlet_v(const dealii::Point<2> & p,int BC_ID)
    {
        for (int i = 0; i<BC_v.size();i++)
        {
            if (BC_ID == BC_v[i])
                return BCval_v[i];
        }
        return 0;
    }
    
    double Dirichlet_p(const dealii::Point<2> & p,int BC_ID)
    {
        for (int i = 0; i<BC_p.size();i++)
        {
            if (BC_ID == BC_p[i])
                return BCval_p[i];
        }
        return 0;
    }
    
    void diff_flux(const Point<2> & p1, const Point<2> & p2,
                 bool &ifconsist, double &valuex, double &valuey)
    {
        ifconsist = false;
        valuex = 0;
        valuey = 0;
    }
    
    void pressure_flux(const Point<2> & p1, const Point<2> & p2,
                   bool &ifconsist, double &valuex, double &valuey)
    {
        ifconsist = false;
        valuex = 0;
        valuey = 0;
    }
    
}


