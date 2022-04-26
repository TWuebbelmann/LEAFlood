# -*- coding: utf-8 -*-
"""
Created on Thu May 20 15:06:37 2021

@author: Wuebbelm
"""

# Import all relevant packages

import cmf
import cmf.geometry 
import numpy as np
import shapefile #PyShp 
import cmf.draw
from shapely.geometry import shape
import pandas as pd 


import os

###

"""
Creates the cmf project
:param subsurface_lateral_connection: Type of lateral subsurface connection, eg. cmf.Darcy
:param surface_lateral_connection:  Type of lateral surface flux connection, eg. cmf.KinematicSurfaceRunoff or None
:param layercount: Number of layers in the subsurface
:return:
"""

def load_meteo(project, meteo_data, begin, timestep, station_name, verbose=True):
    #used by the main function
    """Loads the meteorology from a csv file.
    This is just one example for access your meteo data.
    Thanks to Python, there is no problem to access most type of datasets,
    including relational databases, hdf / netcdf files or what ever else.
    """
    
    if verbose: print ('In: Meteo. Loading Meteorology...')

    # https://philippkraft.github.io/cmf/cmf_tut_test_data.html
    # Create a timeseries for rain - timeseries objects in cmf is a kind of extensible array of
    # numbers, with a begin date, a timestep.
 

    rain = cmf.timeseries(begin=begin, step=timestep)
    
    
    # Create a meteo station
    meteo = project.meteo_stations.add_station(name=station_name, position=(0,0,0))
    # Meteorological timeseries, if you prefer the beginning of the timeseries
    # read in from file, just go ahead, it is only a bit of Python programming
    meteo.Tmax = cmf.timeseries(begin=begin, step=timestep) 
    meteo.Tmin = cmf.timeseries(begin=begin, step=timestep) 
    meteo.rHmean = cmf.timeseries(begin=begin, step=timestep)
    meteo.Windspeed = cmf.timeseries(begin=begin, step=timestep)
    meteo.Sunshine = cmf.timeseries(begin=begin, step=timestep)

    # Load climate data from csv file
    csvfile = open(meteo_data) 
    csvfile.readline()  # Read the headers, and ignore them
    for line in csvfile:
        # split the line in to columns using commas/ semincolons
        #print('Time Step: ',line)
        columns = line.split(';') 
        # Get the values, but ignore the date, we have begin and step
        # of the data file hardcoded
        # If you don't get this line - it is standard Python, I would recommend the official Python.org tutorial
        for timeseries, value in zip([rain, meteo.Tmax, meteo.Tmin, meteo.rHmean, meteo.Windspeed, meteo.Sunshine],
                                     [float(col) for col in columns[1:]]):
            
            # Add the value from the file to the timeseries
            timeseries.add(value)

    meteo.T = (meteo.Tmax + meteo.Tmin) * 0.5  

    # Create a rain gauge station
    project.rainfall_stations.add(station_name, rain, (0, 0, 0))


    # Use the rainfall for each cell in the project
    project.use_nearest_rainfall()
    # Use the meteorological station for each cell of the project
    project.use_nearest_meteo()
    
    if verbose: print ('Reading Meteorology from CSV done')

class Model:
    """
    The 3d model based on irregular polygons from a shapefile
    Set up the Cells with Storages and fluxes
    """
    

    def build_cell(self, c, layercount=1, Ksat = 1, nManning = 0.1, initial_waterlevel=0.,
                   init_sat_depth=0.3, n=0.3, param_b=8., theta_x=0.2, n_decay=0.2):
        #is used by the main function
        # https://philippkraft.github.io/cmf/classcmf_1_1upslope_1_1_cell.html
        """
        Shapes and makes a cell, by adding layers and connections
        :param c: cmf.Cell
        :param layercount: Number of layers in the subsurface (comes with int function)
        :Ksat: saturated conductivity
        :return: the now enhanced cell
        """

        
        d = 0.5 # depth of the soil [m]
        # Linear retention curve. Ksat is defined from the Shapefile for each cell
        # https://philippkraft.github.io/cmf/classcmf_1_1upslope_1_1_linear_retention.html
        
        # _b (Reduction K as a funktion of Ksat) and theta, Potential (Physical -> Hydraulic Head) for droughnes of the soil. Not implemented yet.
        RetentionCurve = cmf.BrooksCoreyRetentionCurve(ksat = Ksat, 
                                                       porosity = n,
                                                       _b=param_b,
                                                       theta_x = theta_x,
                                                       porosity_decay = n_decay)
        
        
        # Add layer of defined thickness and retention curve for each cell
        # https://philippkraft.github.io/cmf/classcmf_1_1upslope_1_1_soil_layer.html
        if layercount > 1:
            for i in range(layercount):
                depth = (i+1) * d/layercount
                c.add_layer(depth,RetentionCurve)
        else:
            c.add_layer(d, RetentionCurve)
        #print('layer & Retention curve created')

        c.saturated_depth = init_sat_depth # set saturated depth: The initial potential of the new layer in m below surface. Default = 10m (=quite dry)
        if c in self.swale_cells:
            c.saturated_depth+=initial_waterlevel
        
        c.surfacewater_as_storage() # define surfacewater as storage
        c.surfacewater.nManning = nManning
        c.surfacewater.puddledepth = 0.01 # minimum depth for runoff
        
        #for c in self.swale_cells:
        #    c.surfacewater.depth = initial_waterlevel #water table depth

              
        
        # Processes within the cell
        c.install_connection(cmf.GreenAmptInfiltration)
        c.install_connection(cmf.RutterInterception) # Interception; https://philippkraft.github.io/cmf/cmf_tut_intercept.html
        cmf.CanopyStorageEvaporation(c.canopy,c.evaporation,c)
   #     c.install_connection(cmf.CanopyOverflow) # does not work yet but should be implemnted

        # Test
        #c.install_connection(cmf.Richards)

        return c
   
        
    

    @staticmethod
    def connect_outlet_cell(outlet_cell: cmf.Cell, subsurface_lateral_connection,
                            surface_lateral_connection):
        o = outlet_cell.surfacewater
        for neighbor, width in outlet_cell.neighbors:
            for l in neighbor.layers:
 #               subsurface_lateral_connection(l, o, width)
 #           try:
                surface_lateral_connection(neighbor.surfacewater, o, width)
 #           except TypeError:
 #               pass

    def get_outflow(self, t=None):
        t = t or self.t
        return sum(o.surfacewater(t) for o in self.outlet_cells) / self.area * 1000
    
  
    #alten Funktion - Liste und Dataframes
    def get_outlet_fluxes(self, t = None):
        """
        Extract the waterbalance at the outlet cell 
        Rows: Cell
        Collumns: Time
        
        returns the dataframes
        
        This can probably be solved in a nicer and more elegant way
        """
        t = t or self.t
        o_id = [] # List with Cell ID
        o_x = [] # List with x-coordinate of cell 
        o_y = [] # List with y-coordinate cell
        o_area = []
            
        #Creating lists for each storage
        list_waterbalance = []
        
        for o in self.outlet_cells:
        # Creating a list for eacht Cell. fill the lists...
            #print(t, '_______',c)
            try:
                # Information of the Cell: ID, Coordinates
                o_id.append(o.Id) 
                o_x.append(o.x)
                o_y.append(o.y)
                o_area.append(o.area)
                
                list_waterbalance.append(o.surfacewater.waterbalance(t)) # waterbalance in m3/day. 
                
            except:
                list_waterbalance.append(np.nan)   
                
                pass
         
        return o_id, o_x, o_y, o_area, list_waterbalance  

    def get_rainfall(self, t=None):
        t = t or self.t
        return sum(c.get_rainfall(t) for c in self.project) / self.area * 1000

    def get_et(self, t=None):
        t = t or self.t
        return sum(c.transpiration(t) + c.evaporation(t) for c in self.project) / self.area * 1000 

    def get_volume(self): # gesamtes Volumen im Gebiet? nicht pro Zelle?
        return sum(sum(s.volume for s in c.storages) for c in self.project) / self.area * 1000 # durch Area*1000? Was ist s?
  
    #used by the run function
    def get_cell_parameters(self, t = None):
        """
        Extract water volume for storages (surface, soil, canopy) for each cell
        Creating dataframes for each storage. 
        Rows: Cell
        Collumns: Time
        
        returns the dataframes
        
        This can probably be solved in a nicer and more elegant way
        """
        
        t = t or self.t
            
        c_id = [] # List with Cell ID
        c_x = [] # List with x-coordinate of cell 
        c_y = [] # List with y-coordinate cell
        c_area = []
            
        #Creating lists for each storage
        list_surfaceDepth = []
        list_soilWater = []
        list_soilWaterCap = []
        list_usedSoilCap = []
        list_interception = []
        # maybe it is usefull to insert the hole water volume of the cell: surface + soil + interception
           
        for c in self.project:
        # Creating a list for eacht Cell. fill the lists...
            try:
                # Information of the Cell: ID, Coordinates
                c_id.append(c.Id) 
                c_x.append(c.x)
                c_y.append(c.y)
                c_area.append(c.area)
                
                # SURFACE Water                        
                list_surfaceDepth.append(c.surfacewater.depth *1000) # returns the water table depth (in m?). converted into mm          
                
                # SOIL Water volume (first layer)  
                soilWaterCont = c.layers.volume[0] / c.area # Returns the volume of water in this storage in m³. Converted in m 
                list_soilWater.append(soilWaterCont *1000) # output in mm
                
                # SOIL water capacity
                l = c.layers[0]
                cap_area = l.get_capacity() / c.area # Returns the capacity of the water storage in m3. Converted in m
                used_cap = soilWaterCont / cap_area
                list_soilWaterCap.append(cap_area) 
                list_usedSoilCap.append(used_cap)
             
                # CANOPY Water 
                list_interception.append((c.canopy.volume / c.area)*1000) #Volume of water stored on the plant [m3]. Converted in m and than to mm
                                
            except:
                list_surfaceDepth.append(np.nan)
                list_soilWater.append(np.nan)
                list_usedSoilCap.append(np.nan)
                list_interception.append(np.nan)
                    
                pass
    
                
        return list_surfaceDepth, list_soilWater, list_usedSoilCap, list_interception

        
    @property
    def t(self): 
        return self.solver.t

    @property
    def area(self):
        return sum(c.area for c in self.project)
    
    # main function. first function, that will open
    def __init__(self, subsurface_lateral_connection,
                       surface_lateral_connection=None,
                       layercount=1, in_shapefile = None, 
                       meteo_data = None, begin = None, 
                       timestep = None, station_name = None, 
                       initial_waterlevel = 0,init_sat_depth=0.3, 
                       n=0.3, param_b=8., theta_x=0.2, n_decay=0.2,
                       scale_nman=1., scale_k=1., no_tree=False,
                       verbose=True):
        """
        Creates the cmf project

        :param subsurface_lateral_connection: Type of lateral subsurface connection, eg. cmf.Darcy
        :param surface_lateral_connection:  Type of lateral surface flux connection, eg. cmf.KinematicSurfaceRunoff or None
        :param layercount: Number of layers in the subsurface
        :return:
        """
        
        self.properties = {'subsurface_lateral_connection': str(subsurface_lateral_connection),
                           'surface_lateral_connection':str(surface_lateral_connection),
                           'layercount':layercount,
                           'shapefile':in_shapefile,
                           'meteo_data':meteo_data,
                           'begin':str(begin),
                           'timestep': str(timestep),
                           'station_name': station_name,
                           'init_w': initial_waterlevel,
                           'init_s': init_sat_depth,
                           'porosity':n,
                           'param_b': param_b,
                           'theta_x': theta_x,
                           'n_decay': n_decay,
                           'scale_nman': scale_nman,
                           'scale_k': scale_k,
                           'no_tree':no_tree}
        
        # create a project
        p = cmf.project()
        if verbose: print('Project created')
        
        # Read the Shapefile with Shapefile from pyshp
        sf = shapefile.Reader(in_shapefile) #Datei wird eingelesen; hier anderes Package für Shapefile, da geos-shapereader nicht mehr läuft
        #valid_sf = make_valid(sf)
        #str(make_valid(sf))
        shp = sf.shapes() # read Geometry
        records = sf.records() # read Attribute
        #feature = sf.shapeRecords()  #store Geometry separately 
        if verbose: print('Read Shapefile: done')
        
        # Create cells
        self.outlet_cells = [] # create list for outlet cells
        self.swale_cells = []
       
        #creating cells based on the shapefile; https://philippkraft.github.io/cmf/classcmf_1_1upslope_1_1_cell.html
        for shp, rec in zip(shp, records):
            # Create a cell for each feature in the shape file. Reads the ID and elevation
            c = cmf.geometry.create_cell(p, shape(shp), rec['Elevation'], rec.oid, with_surfacewater=False) 
            
           # print(c,rec['CLASS'])
            #print (rec.oid) 
            # Defining the LU specific parameters for each Feature
            # If it is an outlet feature, add cell to the list of outletcells
            if rec.CLASS.startswith('Outlet'):
                #print('Outlet Shape: ', rec)
                self.outlet_cells.append(c)
            
           
            else:
                # If it is a normal upload cell, define soil and surface parameters
                #c.surfacewater.nManning = rec['Manning'] #roughness
                Ksat = rec['Ksat']  # saturated conductivity [m/day]
                if Ksat > 5E-5:
                    Ksat*=scale_k
                #c.surfacewater.nManning = rec['Manning']
                nManning = rec['Manning'] * scale_nman
               
                
                if rec.CLASS.startswith('Swale'):
                        if verbose: print(rec['CLASS'], initial_waterlevel)
                        n=0.6
                        self.swale_cells.append(c)
                        

                #Create canopy if available/ needed and set vegetation parameters
                #https://philippkraft.github.io/cmf/cmf_tut_intercept.html
                lai = rec['LAI']
                if no_tree:
                    lai = 0.
                if lai > 0.0:
                    
                    c.add_storage('Canopy','C') #canopy as storage
                    cmf.Rainfall(c.canopy,c,False,True) # RS->canopy, only intercepted rain
                    cmf.Rainfall(c.surfacewater,c,True,False) # RS->surface, only throughfall
                    
                    # https://philippkraft.github.io/cmf/classcmf_1_1upslope_1_1vegetation_1_1_vegetation.html
                    c.vegetation.LAI = lai # rec['LAI']
                    InterceptionCapacity = rec['Int_Cap']
                    c.vegetation.CanopyClosure = rec['Canopy_C']
                    c.vegetation.CanopyCapacityPerLAI = InterceptionCapacity / c.vegetation.LAI

                # if there is no canopy, set canopy closure to 0
                else:
                    #c.surfacewater.nManning = rec['Manning']
                   # c.add_storage('Canopy','C') #canopy as storage
                    #cmf.Rainfall(c.canopy,c,False,True) # RS->canopy, only intercepted rain
                    cmf.Rainfall(c.surfacewater,c,True,False) # RS->surface, only throughfall
                    
                    c.vegetation.CanopyClosure = 0
                    c.vegetation.LAI = rec['LAI']
                    c.vegetation.CanopyCapacityPerLAI = 0
                    
                                    
                # call build cell function from the upper part of the script. Sets the fluxes and soil conditions
                #def build_cell(self, c, layercount=0, Ksat = 1, nManning = 0.1, initial_waterlevel=0.,
                #init_sat_depth=0.3, n=0.3, param_b=8., theta_x=0.2, n_decay=0.2):
                self.build_cell(c, Ksat = Ksat, nManning = nManning, initial_waterlevel=initial_waterlevel,
                                init_sat_depth=init_sat_depth, n=n, param_b=param_b, 
                                theta_x=theta_x, n_decay=n_decay, layercount=layercount) 
                       
        if verbose: print ('Cells created')
        
        # Build topology
        cmf.geometry.mesh_project(p, verbose=verbose)
        if verbose: print ('Geometry created')

        # Connect cells with fluxes | Cell connection with neighbour cells
        #  https://philippkraft.github.io/cmf/namespacecmf.html#a771066e276c39bccafb627cbe1143a6c
        #cmf.connect_cells_with_flux(p, subsurface_lateral_connection) 
        if surface_lateral_connection:
            cmf.connect_cells_with_flux(p, surface_lateral_connection)    
        
        # Connect outlets with neighbor cells
        for o_cell in self.outlet_cells:
            self.connect_outlet_cell(o_cell, subsurface_lateral_connection, surface_lateral_connection)
        
        if verbose: print ('Cells connected')
        
        # Load driver data: Meteorology Data
        if verbose: print('Go to Meteo...')
        load_meteo(p, meteo_data, begin, timestep, station_name, verbose=verbose)
        
        self.project = p
        # I have nothing changed here. But here is the documentation link
        # https://philippkraft.github.io/cmf/cmf_tut_solver.html
        self.solver = cmf.CVodeIntegrator(p, 1e-9)
        self.solver.t = p.meteo_stations[0].T.begin

    
    def run2(self, start=None, end=None, step=cmf.min*5, 
             verbose=True, export_results=True):# Timesteps must be change here. cmf.day
        #called after the main function
        start = start or self.project.meteo_stations[0].T.begin
        end = end or self.project.meteo_stations[0].T.end
        
        if verbose: 
            print('Berechnungszeit: start ', start, ' ; end ', end, '; step ', step)
            
            print('_______________________________________________')
            print('start creating Dataframes per Hydroparameter')
        
        # Creating data frame for the output values from the function 'get_parameter_per_cell'
        list_t = []
        list_i = []
        list_Surfacewater_total = []
        list_Soilwater_total = []
        list_usedSoilwaterCap_total = []
        list_Interception_total = []
        list_waterbalance_total = []
        
        # for each time step
        for i, t in enumerate(self.solver.run(start, end, step)):
            if verbose: print (i, ': ', t)

            list_surfaceDepth, list_soilWater, list_usedSoilCap, list_interception = self.get_cell_parameters(t)
            o_id, o_x, o_y, o_area, list_waterbalance = self.get_outlet_fluxes(t)
            list_t.append(t)
            list_i.append(i)
            
            if i == 0: # first timestep is with cell information (ID and coordinates)
                list_Surfacewater_total.append(list_surfaceDepth) 
                list_Soilwater_total.append(list_soilWater)
                list_usedSoilwaterCap_total.append(list_usedSoilCap)
                list_Interception_total.append(list_interception)
                list_waterbalance_total.append(list_waterbalance)

            else: # append new collumn/ timestep
                list_Surfacewater_total.append(list_surfaceDepth) 
                list_Soilwater_total.append(list_soilWater)
                list_usedSoilwaterCap_total.append(list_usedSoilCap)
                list_Interception_total.append(list_interception)
                list_waterbalance_total.append(list_waterbalance)
        
        df_surfacewater = pd.DataFrame(list_Surfacewater_total)
        dftr_surfacewater = df_surfacewater.transpose()
        dftr_surfacewater.columns = list_i        

        df_SoilWater = pd.DataFrame(list_Soilwater_total)
        dftr_SoilWater = df_SoilWater.transpose()
        dftr_SoilWater.columns = list_i 
        
        df_usedSoilWaterCap = pd.DataFrame(list_usedSoilwaterCap_total)
        dftr_usedSoilWaterCap = df_usedSoilWaterCap.transpose()
        dftr_usedSoilWaterCap.columns = list_i         
        
        df_Interception = pd.DataFrame(list_Interception_total)
        dftr_Interception = df_Interception.transpose()
        dftr_Interception.columns = list_i         
                
        df_waterbalance = pd.DataFrame(list_waterbalance_total,columns=o_id)
        dftr_waterbalance = df_waterbalance.transpose()
        dftr_waterbalance.columns = list_i
        
        
        if export_results:
            if verbose: print('Exporting results...')
            path = r'.'
            folder = 'output'
            
            if not os.path.exists(os.path.join(path, folder)):
                os.makedirs(os.path.join(path, folder))
            
            dftr_surfacewater.to_excel(os.path.join(path, folder, 'SurfaceWater5.xlsx'), float_format='%.5f', na_rep = 'N/A')
            dftr_SoilWater.to_excel(os.path.join(path, folder, 'SoilWaterCont.xlsx'), float_format='%.5f', na_rep = 'N/A')
            dftr_usedSoilWaterCap.to_excel(os.path.join(path, folder, 'usedSoilWaterCap.xlsx'), float_format='%.5f', na_rep = 'N/A')
            dftr_Interception.to_excel(os.path.join(path, folder, 'InterceptWater.xlsx'), float_format='%.5f', na_rep = 'N/A')
            dftr_waterbalance.to_excel(os.path.join(path, folder, 'WaterbalanceOutlets.xlsx'), float_format='%.5f', na_rep = 'N/A')
            if verbose: print('Done.')
        return dftr_waterbalance, dftr_surfacewater, dftr_SoilWater, dftr_Interception # returns time series for each outlet
    
    def get_time_dim(self):
        return self.project.rainfall_stations[0].data.to_pandas().index

if __name__ == '__main__':    

    # input files        
    in_shapefile = r'GIS/sample.shp'
    meteo_data = r'meteo/meteo.csv' 
    
    station_name = 'Station'
    begin = cmf.Time(4,9,2011,13,0)
    
    # time step
    timestep = cmf.min*60
    
    # initial condition (saturated depth)
    init_sat_depth = 0.3
    
    
    #Call Model
    print('Set Up Model...')
    #call the main function
    m = Model(cmf.Darcy, cmf.KinematicSurfaceRunoff, 1, in_shapefile, meteo_data, begin, timestep, station_name, init_sat_depth=init_sat_depth, verbose=False)
    print('Run Model...')
    
    #run the modell and creating output data
    results = m.run2(step=timestep,verbose=False, export_results=True)
    # reformat results for element #2 (outlet, see shape) in m^3 / day
    outflow = pd.Series(results[0].T[2].values, 
                        index=m.get_time_dim())
    print('end')

    # plot outflow
    outflow.plot()

