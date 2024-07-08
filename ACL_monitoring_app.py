

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import webbrowser
import plotly.io as pio
import os
from scipy.signal import find_peaks
from scipy import integrate
from nptdms import TdmsFile
from openpyxl import load_workbook
from datetime import datetime
from dateutil.relativedelta import relativedelta
import webbrowser






def CMJ_analysis(df):
    #New Column Sum of vertical Force
    df["total_fz"]=df["right_fz"]+df["left_fz"]
    
    # Step 1: Calculate Body weight
    #Find period of quiet standing
    segment_length = 1000
    min_variance = float('inf')
    min_variance_segment = None
   
    #Identify a period of quiet standing during the first 2 seconds.
    for i in range(0, 2000):
        segment = df['total_fz'].iloc[i:i+segment_length]
        segment_variance = segment.var()
        if segment_variance < min_variance:
            min_variance = segment_variance
            min_variance_segment = segment
    bodyweight = min_variance_segment.mean()
    SD5 = min_variance_segment.std() *5
    bodyweight_kg = bodyweight /9.81
    
    #STEP 2: Detect Peaks
    #Detect Jumps using instances of flights
    df["normal_fz"] = (df["total_fz"]-bodyweight)*-1  #Create column for total F minus BW, inverted
    flights, _ = find_peaks(df["normal_fz"], prominence=500, distance=5000)   #detect peaks
       
    #Identify start and end from detected peaks
    CMJ_start_list = []
    CMJ_end_list = []
    for i in flights:
        CMJ_start = i - 1900
        CMJ_start_list.append(CMJ_start)
        CMJ_end = i + 1500
        CMJ_end_list.append(CMJ_end)
            
    #Create a dataframe for each jump with label CMJ_i
    jumps = {}
    for i, (start, end) in enumerate(zip(CMJ_start_list, CMJ_end_list), 1):
        # Create a new dataframe for each iteration
        df_name = f'Countermovement Jump {i}'
        jumps[df_name] = df[["total_fz", "left_fz", "right_fz"]][start:end]
    for df_name, df_data in jumps.items():
        # Create a new column called 'time' where values are derived from the index divided by 1000
        df_data['time'] = (df_data.index - df_data.index[0]) / 1000
    
    CMJ_df_analyzed = pd.DataFrame()
    
    for jump_number, jump_df in jumps.items():
        print("Jump Number:", jump_number)
        #Identify Start
        CMJ_start1 = jump_df.index[jump_df["total_fz"] > bodyweight + 25].min()  #(Perez-Castilla et al (2019) found that 5SD of system weight was most relaible threshold while including the most of the force-time signal)
        CMJ_start2 = jump_df.index[jump_df["total_fz"] < bodyweight - 25].min()  
        if CMJ_start1 < CMJ_start2:
            start = CMJ_start1
        else:
            start = CMJ_start2
            
        #Identify takeoff
        takeoff = jump_df.index[jump_df["total_fz"] < 30].min() 
        
        #Identify Braking and Propulsive Phases
        jump_df['acceleration'] = (jump_df["total_fz"] - bodyweight )/ bodyweight_kg
        jump_time = jump_df.loc[start:takeoff, "time"]
        jump_df["velocity"] = integrate.cumtrapz(jump_df["acceleration"], jump_df["time"], initial = 0)
         
        braking_start = (jump_df['velocity'].loc[:takeoff]).idxmin()
        subset_df2 = jump_df['velocity'].loc[braking_start:takeoff]
        braking_end = (subset_df2>0).idxmax()
         
        #Flight Height
        jump_force = jump_df.loc[start:takeoff, "total_fz"]
        jump_time = jump_df.loc[start:takeoff, "time"]
        GRF_impulse = np.trapz(jump_force, x=jump_time)
        BW_force = np.full_like(jump_force, bodyweight)
        BW_impulse = np.trapz(BW_force, x=jump_time)
        
        #Calculate takeoff velocity
        takeoff_velocity = (GRF_impulse - BW_impulse) /bodyweight_kg
        
        #Calculate Flight height
        flight_height = (takeoff_velocity**2)/(9.81*2)
        print('CMJ JH:', flight_height)
         
        #Landing Force    
        R_landing_F = jump_df.loc[takeoff: , "right_fz"].max()
        L_landing_F = jump_df.loc[takeoff: , "left_fz"].max()
        
        #Landing Force Asymmetry
        landing_F_AI = ((L_landing_F - R_landing_F) / R_landing_F * 100) if R_landing_F > L_landing_F else ((L_landing_F - R_landing_F) / L_landing_F * 100)
        
        #Impulse Braking Phase
        braking_duration = jump_df.loc[braking_start:braking_end, "time"].values
        R_braking_force = jump_df.loc[braking_start:braking_end, "right_fz"].values
        R_braking_impulse = np.trapz(R_braking_force, x=braking_duration)
        L_braking_force = jump_df.loc[braking_start:braking_end, "left_fz"].values
        L_braking_impulse= np.trapz(L_braking_force, x=braking_duration)
        
        #Braking Duration
        braking_duration = jump_df.at[braking_end, 'time'] - jump_df.at[braking_start, 'time']

        #Impulse Braking Asymmetry
        braking_impulse_AI = ((L_braking_impulse - R_braking_impulse) / R_braking_impulse * 100) if R_braking_impulse > L_braking_impulse else ((L_braking_impulse - R_braking_impulse) / L_braking_impulse * 100)
        #print("Right braking impulse:", R_braking_impulse, "Left Braking impulse", L_braking_impulse)
        
        #Impulse Propulsive Phase
        propulsive_duration = jump_df.loc[braking_end:takeoff, "time"].values
        R_propulsive_force = jump_df.loc[braking_end:takeoff, "right_fz"].values
        R_propulsive_impulse= np.trapz(R_propulsive_force, x=propulsive_duration)
        L_propulsive_force = jump_df.loc[braking_end:takeoff, "left_fz"].values
        L_propulsive_impulse = np.trapz(L_propulsive_force, x=propulsive_duration)
        
        #Propulsive Duration
        propulsive_duration = jump_df.at[takeoff, 'time'] - jump_df.at[braking_end, 'time']

                #Impulse Propulsive Asymmetry
        propulsive_impulse_AI = ((L_propulsive_impulse - R_propulsive_impulse) / R_propulsive_impulse * 100) if R_propulsive_impulse > L_propulsive_impulse else ((L_propulsive_impulse - R_propulsive_impulse) / L_propulsive_impulse * 100)
         
        #Plot Left vs Right
        fig = go.Figure()
        # Add traces for left and right forces
        fig.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["left_fz"], mode='lines', name='Left', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["right_fz"], mode='lines', name='Right', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["total_fz"], mode='lines', name='Total', line=dict(color='green')))
        # Add vertical lines for takeoff and start
        fig.add_vline(x=jump_df.loc[takeoff, 'time'], line=dict(color='black', dash='dash'), name='Takeoff', annotation_position="top", annotation_text="Takeoff")
        fig.add_vline(x=jump_df.loc[start, 'time'], line=dict(color='black', dash='dash'), name='Start',  annotation_position="top", annotation_text="Start")
        # Add shaded regions for braking and propulsive phases
        fig.add_vrect(x0=jump_df.loc[braking_start, 'time'], x1=jump_df.loc[braking_end, 'time'], fillcolor="grey", opacity=0.1, line_width=0, annotation_text="Braking Phase", annotation_position="top left", annotation_textangle=-90)
        fig.add_vrect(x0=jump_df.loc[braking_end, 'time'], x1=jump_df.loc[takeoff, 'time'], fillcolor="darkgrey", opacity=0.1, line_width=0, annotation_text="Propulsive Phase", annotation_position="top left", annotation_textangle=-90)

        # Update layout
        fig.update_layout(
            title=f'Jump {jump_number}',
            xaxis_title='Seconds',
            yaxis_title='Force (N)',
            legend=dict(x=0, y=1),
                                    )
        #pio.write_html(fig, file=f'{jump_number}.html', auto_open=True)
        st.plotly_chart(fig)
        
        df_row = pd.DataFrame({ 
            'Bodyweight (kg)': [bodyweight_kg],
            'CMJ Jump Number': [jump_number],
            'CMJ Jump Height (m)': [flight_height],
            'CMJ Right Braking Impulse (Ns)': [R_braking_impulse],
            'CMJ Left Braking Impulse (Ns)': [L_braking_impulse],
            'CMJ Right Braking Relative Impulse (Ns/kg)': [R_braking_impulse/bodyweight_kg],
            'CMJ Left Braking Realtive Impulse (Ns/kg)': [L_braking_impulse/bodyweight_kg],
            'CMJ Braking Impulse Asymmetry Index (%)': [braking_impulse_AI],
            'CMJ Braking Phase Duration (s)': [braking_duration],
            'CMJ Right Propulsive Impulse (Ns)': [R_propulsive_impulse],
            'CMJ Left Propulsive Impulse (Ns)': [L_propulsive_impulse],
            'CMJ Right Propulsive Relative Impulse (Ns/kg)': [R_propulsive_impulse/bodyweight_kg],
            'CMJ Left Propulsive Relative Impulse (Ns/kg)': [L_propulsive_impulse/bodyweight_kg],
            'CMJ Propulsive Impulse Asymmetry Index (%)': [propulsive_impulse_AI],
            'CMJ Propulsive Phase Duration (s)' : [propulsive_duration],
            'CMJ Right Peak Landing Force (N)': [R_landing_F],
            'CMJ Left Peak Landing Force (N)': [L_landing_F],
            'CMJ Right Relative Peak Landing Force (N/kg)': [R_landing_F/bodyweight_kg],
            'CMJ Left Relative Peak Landing Force (N/kg)': [L_landing_F/bodyweight_kg],
            'CMJ Landing Peak Force Asymmetry (%)': [landing_F_AI]  })
         
        CMJ_df_analyzed = pd.concat([CMJ_df_analyzed, df_row], ignore_index=True)    
    
    return CMJ_df_analyzed
    
def DJ_analysis(df):
    #New Column Sum of vertical Force
    df["total_fz"]=df["right_fz"]+df["left_fz"]
        

    # Detect peaks
    flights, _ = find_peaks(df["total_fz"],
                            height=2000,
                            threshold=None,
                            distance=5000,
                            prominence=None,
                            width=None,
                            plateau_size=None)
  
    #Identify Start and End Indexes
    DJ_start_list = []
    DJ_end_list = []
    for i in flights:
        DJ_start = i - 1500
        DJ_start_list.append(DJ_start)
        DJ_end = i + 2400
        DJ_end_list.append(DJ_end)
    #Ensure start and end indexes are correct

    jumps = {}
    for i, (start, end) in enumerate(zip(DJ_start_list, DJ_end_list), 1):
        df_name = f'Drop Jump {i}'
        jumps[df_name] = df.loc[start:end, ["total_fz", "left_fz", "right_fz"]].copy()

    #Create empty list for each variable. Will Append these below
    DJ_bodyweight_list = []
    DJ_df_analyzed= pd.DataFrame()
        
    for jump_number, jump_df in jumps.items():
        st.write("Jump Number:", jump_number)
        
        jump_df['time'] = (jump_df.index - jump_df.index[0]) / 1000
        #Identify key instants
        ground_contact = (jump_df["total_fz"]>30).idxmax()
        jump_df = jump_df.loc[ground_contact:].copy()
        takeoff = (jump_df["total_fz"].loc[ground_contact+100:]<30).idxmax()
        landing = (jump_df["total_fz"].loc[takeoff+100:]>30).idxmax()        
       
        # Calculate BW for each jump
        bodyweight = jump_df['total_fz'].tail(500).mean()
        bodyweight_kg = bodyweight / 9.81
        DJ_bodyweight_list.append(bodyweight_kg)
       
        #Calculate velocity by integrating acceleration then offsetting by initial velocity. 
        jump_df["acceleration"] = (jump_df["total_fz"] - bodyweight) / bodyweight_kg
        jump_df["velocity"] = integrate.cumtrapz(jump_df["acceleration"], jump_df["time"], initial=0)
        jump_df["velocity"] = jump_df["velocity"] - 2.426
        
        velocity_offset = jump_df['velocity'].tail(500).mean()
                
        jump_df["velocity_post_offset"] = jump_df["velocity"] - velocity_offset
        
        #Calculate jump height using impulse momentum
        jump_force = jump_df.loc[ground_contact:takeoff, "total_fz"]
        jump_time = jump_df.loc[ground_contact:takeoff, "time"]
        GRF_impulse = np.trapz(jump_force, x=jump_time)
        BW_force = np.full_like(jump_force, bodyweight)
        BW_impulse = np.trapz(BW_force, x=jump_time)
     
        #Calculate takeoff velocity
        takeoff_velocity = ((GRF_impulse - BW_impulse) - (bodyweight_kg*(2.426-velocity_offset))) /bodyweight_kg
        takeoff_velocity2 = ((GRF_impulse - BW_impulse) - (bodyweight_kg*(2.426))) /bodyweight_kg

        #Calculate Flight height
        flight_time = jump_df.at[landing, 'time'] - jump_df.at[takeoff, 'time']
        flight_height_flight_time_method = (9.81*flight_time**2)/8
        flight_height_impulse_method = (takeoff_velocity**2)/(9.81*2)
        flight_height_impulse_method2 = (takeoff_velocity2**2)/(9.81*2)
        fall_height_postoffset = ((-2.462-velocity_offset)**2)/(19.62)

        #Calculate RSI and append to list
        ground_contact_time = jump_df.at[takeoff, 'time'] - jump_df.at[ground_contact, 'time']
        RSI = flight_height_impulse_method / ground_contact_time
        
        
        #Calculate End of Braking  / Start of Propulsive         
        subset_df = jump_df['velocity_post_offset'].loc[ground_contact:takeoff]
        braking_end = (subset_df>0).idxmax()
        
        #Calculate Left and Right Braking Impulset
        braking_duration = jump_df.loc[ground_contact:braking_end, "time"].values
        R_braking_force = jump_df.loc[ground_contact:braking_end, "right_fz"].values
        R_braking_impulse = np.trapz(R_braking_force, x=braking_duration)
        L_braking_force = jump_df.loc[ground_contact:braking_end, "left_fz"].values
        L_braking_impulse = np.trapz(L_braking_force, x=braking_duration)
        
        #Impulse Braking Asymmetry
        braking_impulse_AI = ((L_braking_impulse - R_braking_impulse) / R_braking_impulse * 100) if R_braking_impulse > L_braking_impulse else ((L_braking_impulse - R_braking_impulse) / L_braking_impulse * 100)
        
        #Impulse Propulsive Phase
        propulsive_duration = jump_df.loc[braking_end:takeoff, "time"].values
        R_propulsive_force = jump_df.loc[braking_end:takeoff, "right_fz"].values
        R_propulsive_impulse = np.trapz(R_propulsive_force, x=propulsive_duration)
        L_propulsive_force = jump_df.loc[braking_end:takeoff, "left_fz"].values
        L_propulsive_impulse = np.trapz(L_propulsive_force, x=propulsive_duration)
        
        #Impulse Propulsive Asymmetry
        propulsive_impulse_AI = ((L_propulsive_impulse - R_propulsive_impulse) / R_propulsive_impulse * 100) if R_propulsive_impulse > L_propulsive_impulse else ((L_propulsive_impulse - R_propulsive_impulse) / L_propulsive_impulse * 100)
        
        
        #Calculate Landing Force for each leg and append to list
        L_landing_force = (jump_df["left_fz"].loc[landing:]).max()
        R_landing_force = (jump_df["right_fz"].loc[landing:]).max()
        landing_F_AI = ((L_landing_force - R_landing_force) / R_landing_force * 100) if R_landing_force > L_landing_force else ((L_landing_force - R_landing_force) / L_landing_force * 100)
        
        #Plot Left vs Right
        fig = go.Figure()
        # Add traces for left and right forces
        fig.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["left_fz"], mode='lines', name='Left', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["right_fz"], mode='lines', name='Right', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["total_fz"], mode='lines', name='Total', line=dict(color='green')))
        # Add vertical lines for takeoff
        fig.add_vline(x=jump_df.loc[takeoff, 'time'], line=dict(color='black', dash='dash'), name='Takeoff', annotation_position="top", annotation_text="Takeoff")
 
        # Update layout
        fig.update_layout(
            title=f'Jump {jump_number}',
            xaxis_title='Seconds',
            yaxis_title='Force (N)',
            legend=dict(x=0, y=1),
                                    )
        st.plotly_chart(fig)

        fig1 = go.Figure()
        # Add traces for left and right forces
        fig1.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["velocity"], mode='lines', name='Pre Offset', line=dict(color='dark grey')))
        fig1.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["velocity_post_offset"], mode='lines', name='Post Offset', line=dict(color='grey')))
        # Add vertical lines for takeoff and start
        fig1.add_hline(y=0, line=dict(color='black', dash='dash'))
      
        # Update layout
        fig1.update_layout(
            title=f'Jump {jump_number} Pre and Post Velocity Offset',
            xaxis_title='Seconds',
            yaxis_title='Velocity (m/s)',
            legend=dict(x=0, y=1),
                                        )
        st.plotly_chart(fig1)

        st.write('Jump Height (Flight Time Method) (m):', flight_height_flight_time_method)
        st.write('Jump Hieght (Impulse-Momentum Method w/ Offset) (m):',flight_height_impulse_method)
        st.write('Jump Height (Impulse-Momentum Method w/out Offset) (m)',flight_height_impulse_method2)
        st.write("Actual Fall Height (m):", fall_height_postoffset)


        df_row = pd.DataFrame({
         'DJ Jump Number': [jump_number],
         'DJ Jump Height (m)': [flight_height_impulse_method],
        'GCT (s)': [ground_contact_time],
        'RSI': [RSI],
        'Actual Fall Height (m)': [fall_height_postoffset],
        'DJ Right Braking Impulse (Ns)': [R_braking_impulse],
        'DJ Left Breaking Impulse (Ns)': [L_braking_impulse],
        'DJ Braking Impulse Asymmetry (%)' : [ braking_impulse_AI],
        'DJ Right Propulsive Impulse (Ns)': [R_propulsive_impulse],
        'DJ Left Propulsive Impulse (Ns)': [L_propulsive_impulse],
        'DJ Propulsive Impulse Asymmetry (%)' : [propulsive_impulse_AI],
        'DJ Right Landing Force (N)': [R_landing_force],
        'DJ Left Landing Force (N)': [L_landing_force],
        'DJ Landing Force Asymmetry Index (%)': [landing_F_AI]  })      

        DJ_df_analyzed = pd.concat([DJ_df_analyzed, df_row], ignore_index=True)
    
    return(DJ_df_analyzed)


def SJ_analysis(df):
    #New Column Sum of vertical Force
    df["total_fz"]=df["right_fz"]+df["left_fz"]

    # Step 1: Calculate Body weight
    #Find period of quiet standing
    segment_length = 1000
    min_variance = float('inf')
    min_variance_segment = None
     #Identify a period of quiet standing during the first 2 seconds.
    for i in range(0, 2000):
        segment = df['total_fz'].iloc[i:i+segment_length]
        segment_variance = segment.var()
        if segment_variance < min_variance:
             min_variance = segment_variance
             min_variance_segment = segment
    bodyweight = min_variance_segment.mean()
    bodyweight_kg = bodyweight /9.81

    #STEP 2: Detect Peaks
    #Detect Jumps using instances of flights
    df["normal_fz"] = (df["total_fz"]-bodyweight)*-1  #Create column for total F minus BW, inverted
    flights, _ = find_peaks(df["normal_fz"], prominence=500, distance=5000)   #detect peaks

    #Identify start and end window from detected peaks
    SJ_start_list = []
    SJ_end_list = []
    for i in flights:
         SJ_start = i - 1000
         SJ_start_list.append(SJ_start)
         SJ_end = i + 1500
         SJ_end_list.append(SJ_end)
     
    #Create a dataframe for each jump with label CMJ_i
    jumps = {}
    for i, (start, end) in enumerate(zip(SJ_start_list, SJ_end_list), 1):
         # Create a new dataframe for each iteration
         df_name = f'Squat Jump {i}'
         jumps[df_name] = df[["total_fz", "left_fz", "right_fz"]][start:end]

    # Create a new column for each called 'time' where values are derived from the index divided by 1000
    for df_name, df_data in jumps.items():
         df_data['time'] = (df_data.index - df_data.index[0]) / 1000
     

    SJ_df_analyzed= pd.DataFrame()

    for jump_number, jump_df in jumps.items():
        #Identify Start. Threshold is 25 N. Other options too variable in SJ.       
        SJ_start1 = jump_df.index[jump_df["total_fz"] > bodyweight+25].min()  #(Perez-Castilla et al (2019) found that 5SD of system weight was most relaible threshold while including the most of the force-time signal)
        SJ_start2 = jump_df.index[jump_df["total_fz"] < bodyweight-25].min()  
        if SJ_start1 < SJ_start2:
            start = SJ_start1
        else:
            start = SJ_start2
            
        #Check for micro countermovement. If present, disregard this jump    
        threshold = bodyweight - (bodyweight * 0.025)
        if (jump_df.loc[start:start-500, 'total_fz'] < threshold).any():
            print(f"Microcountermovement identified on jump # {jump_number} ")
            continue
        print("No MicroCM detected")                    
        #Identify takeoff
        takeoff = jump_df.index[jump_df["total_fz"] < 30].min() 
        
        #Flight Height
        jump_force = jump_df.loc[start:takeoff, "total_fz"]
        jump_time = jump_df.loc[start:takeoff, "time"]
        GRF_impulse = np.trapz(jump_force, x=jump_time)
        BW_force = np.full_like(jump_force, bodyweight)
        BW_impulse = np.trapz(BW_force, x=jump_time)
        
        #Calculate takeoff velocity
        takeoff_velocity = (GRF_impulse - BW_impulse) /bodyweight_kg
        
        #Calculate Flight height
        flight_height = (takeoff_velocity**2)/(9.81*2)
        #print("Jump height: ", flight_height)
         
        #Landing Force    
        R_landing_F = jump_df.loc[takeoff: , "right_fz"].max()
        L_landing_F = jump_df.loc[takeoff: , "left_fz"].max()
         
        #Landing Force Asymmetry
        landing_F_AI = ((L_landing_F - R_landing_F) / R_landing_F * 100) if R_landing_F > L_landing_F else ((L_landing_F - R_landing_F) / L_landing_F * 100)
        
        #Identify Peak Force for Early and Late Phase analysis
        peak_force = (jump_df["total_fz"].loc[start:takeoff]).idxmax()

        
        #Impulse Early Phase
        early_duration = jump_df.loc[start:peak_force, "time"].values
        R_early_force = jump_df.loc[start:peak_force, "right_fz"].values
        R_early_impulse = np.trapz(R_early_force, x=early_duration)
        L_early_force = jump_df.loc[start:peak_force, "left_fz"].values
        L_early_impulse = np.trapz(L_early_force, x=early_duration)
        
        #Impulse Braking Asymmetry
        early_impulse_AI = ((L_early_impulse - R_early_impulse) / R_early_impulse * 100) if R_early_impulse > L_early_impulse else ((L_early_impulse - R_early_impulse) / L_early_impulse * 100)
        
        #Impulse Propulsive Phase
        late_duration = jump_df.loc[peak_force:takeoff, "time"].values
        R_late_force = jump_df.loc[peak_force:takeoff, "right_fz"].values
        R_late_impulse = np.trapz(R_late_force, x=late_duration)
        L_late_force = jump_df.loc[peak_force:takeoff, "left_fz"].values
        L_late_impulse = np.trapz(L_late_force, x=late_duration)
        
        #Impulse Propulsive Asymmetry
        late_impulse_AI = ((L_late_impulse - R_late_impulse) / R_late_impulse * 100) if R_late_impulse > L_late_impulse else ((L_late_impulse - R_late_impulse) / L_late_impulse * 100)
         
        #Plot Left vs Right
        fig = go.Figure()
        # Add traces for left and right forces
        fig.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["left_fz"], mode='lines', name='Left', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["right_fz"], mode='lines', name='Right', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=jump_df['time'], y=jump_df["total_fz"], mode='lines', name='Total', line=dict(color='green')))
        # Add vertical lines for takeoff and start
        fig.add_vline(x=jump_df.loc[takeoff, 'time'], line=dict(color='black', dash='dash'), name='Takeoff', annotation_position="top", annotation_text="Takeoff")
        fig.add_vline(x=jump_df.loc[start, 'time'], line=dict(color='black', dash='dash'), name='Start',  annotation_position="top", annotation_text="Start")
        # Add shaded regions for early and late phases
        fig.add_vrect(x0=jump_df.loc[start, 'time'], x1=jump_df.loc[peak_force, 'time'], fillcolor="grey", opacity=0.1, line_width=0, annotation_text="Early Phase", annotation_position="top left", annotation_textangle=-90)
        fig.add_vrect(x0=jump_df.loc[peak_force, 'time'], x1=jump_df.loc[takeoff, 'time'], fillcolor="darkgrey", opacity=0.1, line_width=0, annotation_text="Late Phase", annotation_position="top left", annotation_textangle=-90)

        # Update layout
        fig.update_layout(
            title=f'Jump {jump_number}',
            xaxis_title='Seconds',
            yaxis_title='Force (N)',
            legend=dict(x=0, y=1),
                                    )
        st.plotly_chart(fig)
        
        df_row = pd.DataFrame({ 
            'SJ Jump Number': [jump_number],
            'SJ Jump Height (m)': [flight_height],
         'Left Early Phase Impulse (Ns)': [L_early_impulse],
         'Right Early Phase Impulse (Ns)': [R_early_impulse],
         'Left Early Phase Relative Impulse (Ns/kg)': [L_early_impulse/bodyweight_kg],
         'Right Early Phase Relative Impulse (Ns/kg)': [R_early_impulse/bodyweight_kg],
         'Early Phase Impulse Asymmetry Index (%)': [early_impulse_AI],
         'Left Late Phase Impulse (Ns)': [L_late_impulse],
         'Right Late Phase Relative Impulse (Ns/kg)': [R_late_impulse/bodyweight_kg],
         'Left Late Phase Relative Impulse (Ns/kg)': [L_late_impulse/bodyweight_kg],
         'Right Late Phase Impulse (Ns)': [R_late_impulse],
         'Late Phase Impulse Asymmetry Index (%)': [late_impulse_AI],         
         'SJ Right Landing Force (N)': [R_landing_F],
         'SJ Left Landing Force (N)': [L_landing_F],
         'SJ Right Relative Landing Force (N/kg)': [R_landing_F/bodyweight_kg],
         'SJ Left Relative Landing Force (N/kg)': [L_landing_F/bodyweight_kg],
         'SJ Landing Force Asymmetry Index (%)': [landing_F_AI]  }) 

        SJ_df_analyzed = pd.concat([SJ_df_analyzed, df_row], ignore_index=True)
 
    return SJ_df_analyzed


def readforceplatetdms(file):
    read_tdms = TdmsFile.read(file) 
    temp_tdms_frame = read_tdms.as_dataframe()
    df = pd.DataFrame(temp_tdms_frame)
    column_names =['left_fx', 'left_fy', 'left_fz', 'left_mx', 'left_my', 'left_mz', 'right_fx', 'right_fy', 'right_fz', 'right_mx', 'right_my', 'right_mz']
    df = df.set_axis(column_names, axis="columns")
    return df


st.session_state['stage'] = 0
@st.experimental_dialog("Upload Jumps:", width='large')
def modal():
    def set_state(i):
        st.session_state['stage'] = i
    
    if st.session_state.stage >= 0:
        CMJ_file = st.file_uploader(label="Select a Countermovement Jump", key=4)
        SJ_file = st.file_uploader(label="Select a Squat Jump", key=5)
        DJ_file = st.file_uploader(label="Select a Drop Jump", key=6)
       
        input_athlete_name = None
        test_date = None
        jump_type_correct = None 
       
        if CMJ_file and SJ_file and DJ_file:
             #Get Athlete Name
             CMJ_file_name = CMJ_file.name
             CMJ_athlete_name = CMJ_file_name.split('_')[1]
             SJ_file_name = SJ_file.name
             SJ_athlete_name = SJ_file_name.split('_')[1]
             DJ_file_name = DJ_file.name
             DJ_athlete_name = DJ_file_name.split('_')[1]
             if DJ_athlete_name == CMJ_athlete_name and CMJ_athlete_name == SJ_athlete_name:
                 input_athlete_name = CMJ_athlete_name
             else:
                 st.header('Error: Please ensure all tests are from the same athlete')
                 input_athlete_name = None
             
             #And Test Date
             CMJ_date_component = CMJ_file_name.split('_')[0]
             CMJ_test_date = datetime.strptime(CMJ_date_component, "%d%m%Y").date()
             
             SJ_date_component = SJ_file_name.split('_')[0]
             
             DJ_date_component = DJ_file_name.split('_')[0]
             
             if CMJ_date_component == DJ_date_component and CMJ_date_component == SJ_date_component:
                 test_date = CMJ_test_date
             else:
                 st.header('Error: Please ensure all tests are from the same date')
                 test_date = None
                 
             #And Test Type
             CMJ_jump_type = CMJ_file_name.split('_')[2]
             SJ_jump_type = SJ_file_name.split('_')[2]
             DJ_jump_type = DJ_file_name.split('_')[2]
             if CMJ_jump_type == 'CMJ' and SJ_jump_type == 'SJ' and DJ_jump_type == 'DJ':
                 jump_type_correct = True
             else: 
                 st.header('Error: Please ensure jumps are uploaded to correct file uploader')
                 jump_type_correct = None

            
            
            
             if input_athlete_name and test_date and jump_type_correct:
                 st.button('Analyze', on_click=set_state, args=[1], key='analyze_button')
             else:
                 st.button('Analyze', disabled=True, key='disabled_analyze_button')
                 
        else: 
            st.button('Please Upload Jumps', disabled=True, key='please_upload_jumps_button')
    
    
    if st.session_state.stage == 1:
        
        st.button('Save', on_click=set_state, args=[2], key='save_button')

        if CMJ_file:
            st.subheader('Countermovement Jumps', divider='gray')
            CMJ_df = readforceplatetdms(CMJ_file)
            st.session_state['CMJ_df_analyzed'] = CMJ_analysis(CMJ_df)
            st.dataframe(st.session_state['CMJ_df_analyzed'])

                    
        if SJ_file:
            st.subheader('Countermovement Jumps', divider='gray')
            SJ_df = readforceplatetdms(SJ_file)
            st.session_state['SJ_df_analyzed'] = SJ_analysis(SJ_df)
            st.dataframe(st.session_state['SJ_df_analyzed'])

        
        if DJ_file:
            st.subheader('Squat Jumps', divider='gray')
            DJ_df = readforceplatetdms(DJ_file)
            st.session_state['DJ_df_analyzed'] = DJ_analysis(DJ_df)
            st.dataframe(st.session_state['DJ_df_analyzed'])

    if st.session_state.stage == 2:
        found=False
               
        if input_athlete_name and test_date:
            st.subheader(f"Athlete Name: {input_athlete_name}")
            database_spreadsheet = 'C:/Users/kiera/OneDrive/Documents/School/Strength and Conditioning/Vikes Research/Forceplate Analysis/ACLR AthleteInfo Database.xlsx'
            book = load_workbook(database_spreadsheet)
            sheet = book.active
            for cell in sheet['A']:
                if cell.value == input_athlete_name:
                        found = True
                        sport = sheet[f'F{cell.row}'].value
                        ACLR_limb = sheet[f'G{cell.row}'].value
                        ACLR_date = sheet[f'C{cell.row}'].value
                        graft_type = sheet[f'H{cell.row}'].value
                        st.write(f"{input_athlete_name} already has information in database:")
                        st.write('Sport:', sport)
                        st.write('Graft Type:', graft_type)
                        st.write('Surgical Limb:', ACLR_limb)
                        break
            if not found:
                st.write("Athlete not yet in database. Please Enter the Following Information")
                    
                sport = st.selectbox('Sport', ("WRUG", "WFH"), placeholder="Sport...", key='sport_for_save_unique', index=None)
                ACLR_limb = st.selectbox('Limb', ("Right", "Left"), placeholder="ACLR Limb...", key='limb_for_save_unique', index=None)
                graft_type = st.selectbox('Graft Type', ("QT", "HT"), placeholder="Graft Type...", key='graft_type_for_save_unqiue', index=None)
                ACLR_date = st.date_input('Surgery Date', format='MM/DD/YYYY', key = 'ACLR_date_for_save')        
                    
                # Calculate the difference using relativedelta
            test_date = datetime.strptime(CMJ_date_component, "%d%m%Y").date()
            test_date = pd.to_datetime(test_date)  # Ensure test_date is a datetime object
            ACLR_date = pd.to_datetime(ACLR_date)  # Ensure ACLR_date is a datetime object
            
            difference = relativedelta(test_date, ACLR_date)
            time_since_surgery = round(difference.years * 12 + difference.months + difference.days / 30.44, 1 )
                                                
            #Create Dataframe with All the info
            info_df = pd.DataFrame({
                    'Athlete Name' : [input_athlete_name, input_athlete_name, input_athlete_name],
                    'Testdate (mm/dd/yyyy)' : [test_date, test_date, test_date],
                    "ACLR Date (mm/dd/yyyy)" : [ACLR_date, ACLR_date, ACLR_date],
                    "Time Since Surgery (months)" : [time_since_surgery,time_since_surgery,time_since_surgery],
                    'Sport' : [sport,sport, sport],
                    "Injury Side" : [ACLR_limb, ACLR_limb, ACLR_limb],
                    "Graft Type": [graft_type, graft_type, graft_type] }) 
                
            info_df['Testdate (mm/dd/yyyy)'] = info_df['Testdate (mm/dd/yyyy)'].dt.strftime('%m/%d/%Y')
            info_df['ACLR Date (mm/dd/yyyy)'] = info_df['ACLR Date (mm/dd/yyyy)'].dt.strftime('%m/%d/%Y')
                
            if CMJ_file and SJ_file and DJ_file:
                #Export all Analyzed Info to Database
                #Combine above dataframe with jumps to single dataframes and list of dataframes
                existing_df = pd.read_excel(database_spreadsheet)
               
                row_exists = ((existing_df['Athlete Name'] == input_athlete_name) & (existing_df['Testdate (mm/dd/yyyy)'] == test_date.strftime('%m/%d/%Y'))).any()
                
                master_df = pd.concat([info_df, st.session_state['CMJ_df_analyzed'], st.session_state['SJ_df_analyzed'], st.session_state['DJ_df_analyzed']], axis=1)
                 
                st.header('Testing Results')
                st.dataframe(master_df)
            
                if row_exists is True:
                    st.button('Testing data already in ACLR database', disabled=True, key='error_button2')
                elif not (sport and graft_type and ACLR_limb):
                    st.button('Please fill out athlete info', disabled=True, key='error_button1')
                else:
                    if st.button('Save to Database', key='save_to_database'):
                        updated_df = pd.concat([existing_df, master_df], axis=0)
                        updated_df.to_excel(database_spreadsheet, index=False)
                        st.write('Successfully saved. You can close this window.')
                            
        else: 
                st.write('Countermovement Jump, Squat Jump, and Drop Jump Required to Save to Database')
                      
           



    
             
         

#Define plotting fucntion
def barplot_asymmetry_data (df, plot_title, unit, left_value_column, right_value_column, asymmetry_column):
    
   # Filter and aggregate the dataframe
    filtered_df = df[['Time Since Surgery (months)', left_value_column, right_value_column, asymmetry_column]]
    mean_df = filtered_df.groupby('Time Since Surgery (months)').mean().reset_index()
    mean_df = mean_df.round(1)
    std_df = filtered_df.groupby('Time Since Surgery (months)').std().reset_index()
    std_df = std_df.round(1)

    # Create bar traces
    trace_left = go.Bar(
        x=mean_df["Time Since Surgery (months)"],
        y=mean_df[left_value_column],
        name="Left",        
        text=mean_df[left_value_column],
        textposition='inside',
        insidetextanchor='start',
        textangle=0,
        error_y=dict( type='data', array=std_df[left_value_column], visible=True ))
    
    trace_right = go.Bar(
        x=mean_df["Time Since Surgery (months)"],
        y=mean_df[right_value_column],
        name="Right",
        marker=dict(color="#ffe476"),text=mean_df[right_value_column],
        textposition='inside',
        insidetextanchor='start',
        textangle=0,
        error_y=dict( type='data', array=std_df[right_value_column], visible=True ))
    
    trace_asymmetry = go.Scatter(
        x=mean_df["Time Since Surgery (months)"],
        y=mean_df[[left_value_column, right_value_column]].max(axis=1) * 1.18,
        mode='text',
        text=[f"AI: {value}%" for value in mean_df[asymmetry_column]],
        textposition="top center",
        showlegend=False,
            textfont=dict(
            family="Arial",
            size=12,
            color="black"))
    
    # Create layout
    layout = go.Layout(
        title=plot_title,
        xaxis=dict(title="Months Since Surgery", tickvals=mean_df["Time Since Surgery (months)"],),
        yaxis=dict(title=unit),
        barmode='group',
        legend=dict(x=0, y=-0.2, orientation='h'))

    # Create figure
    fig = go.Figure(data=[trace_left, trace_right, trace_asymmetry], layout=layout)
    return fig

def barplot_RSI (df, plot_title, unit, value_column):
    
   # Filter and aggregate the dataframe
    filtered_df = df[['Time Since Surgery (months)', value_column]]
    mean_df = filtered_df.groupby('Time Since Surgery (months)').mean().reset_index()
    mean_df = mean_df.round(2)
    std_df = filtered_df.groupby('Time Since Surgery (months)').std().reset_index()
    std_df = std_df.round(2)
   
    # Create bar traces
    trace = go.Bar(
        x=mean_df["Time Since Surgery (months)"],
        y=mean_df[value_column],
        name="Left",
        text=mean_df[value_column],
        textposition='inside',
        insidetextanchor='start',
        textangle=0,
        error_y=dict( type='data', array=std_df[value_column], visible=True ))
   
    # Create layout
    layout = go.Layout(
        title=plot_title,
        xaxis=dict(title="Months Since Surgery", tickvals=mean_df["Time Since Surgery (months)"],),
        yaxis=dict(title=unit,  range=[0, 1.8]),
        barmode='group',
        legend=dict(x=0, y=-0.2, orientation='h'))

    # Create figure
    fig = go.Figure(data=trace, layout=layout)
    return fig
def barplot_GCT (df, plot_title, unit, value_column):
     
   # Filter and aggregate the dataframe
    filtered_df = df[['Time Since Surgery (months)', value_column]]
    mean_df = filtered_df.groupby('Time Since Surgery (months)').mean().reset_index()
    mean_df = mean_df.round(2)
    std_df = filtered_df.groupby('Time Since Surgery (months)').std().reset_index()
    std_df = std_df.round(2)
    # Create bar traces
    trace = go.Bar(
        x=mean_df["Time Since Surgery (months)"],
        y=mean_df[value_column],
        name="Left",
        text=mean_df[value_column],
        textposition='inside',
        insidetextanchor='start',
        textangle=0,
        error_y=dict( type='data', array=std_df[value_column], visible=True ))
   
    # Create layout
    layout = go.Layout(
        title=plot_title,
        xaxis=dict(title="Months Since Surgery", tickvals=mean_df["Time Since Surgery (months)"],),
        yaxis=dict(title=unit,  range=[0, 0.6]),
        barmode='group',
        legend=dict(x=0, y=-0.2, orientation='h'))

    # Create figure
    fig = go.Figure(data=trace, layout=layout)
    return fig
def barplot_jump_height (df, plot_title, unit, value_column):
       
   # Filter and aggregate the dataframe
    filtered_df = df[['Time Since Surgery (months)', value_column]]
    mean_df = filtered_df.groupby('Time Since Surgery (months)').mean().reset_index()
    mean_df = mean_df.round(2)
    std_df = filtered_df.groupby('Time Since Surgery (months)').std().reset_index()
    std_df = std_df.round(2)  

    # Create bar traces
    trace = go.Bar(
        x=mean_df["Time Since Surgery (months)"],
        y=mean_df[value_column],
        name="Left",
        text=mean_df[value_column],
        textposition='inside',
        insidetextanchor='start',
        textangle=0,
        error_y=dict( type='data', array=std_df[value_column], visible=True ))
   
    # Create layout
    layout = go.Layout(
        title=plot_title,
        xaxis=dict(title="Months Since Surgery", tickvals=mean_df["Time Since Surgery (months)"],),
        yaxis=dict(title=unit,  range=[0, 0.6]), 
        barmode='group',
        legend=dict(x=0, y=-0.2, orientation='h'))

    # Create figure
    fig = go.Figure(data=trace, layout=layout)
    return fig



url = 'https://raw.githubusercontent.com/kieransphillips/ACL_monitoring/main/ACLR%20AthleteInfo%20Database.xlsx'
df = pd.read_excel(url, sheet_name='Sheet1')

#Streamlit dash
st.set_page_config(layout="wide")


with st.sidebar:
    #Filter to only athlete of interest 
    athlete_names = df['Athlete Name'].unique()
    athlete_name = st.selectbox('Select an Athlete', athlete_names, index=None, placeholder='Select Athlete...')
    filtered_data = df[df["Athlete Name"] == f"{athlete_name}"]
    
    #Select or autofill ACL info
    if athlete_name is not None:
        sport = st.selectbox("Sport", filtered_data["Sport"].unique(), key='sport_filled')
        graft_type = st.selectbox('Graft Type', filtered_data['Graft Type'].unique(), key='graft_type_filled')
        ACLR_limb = st.selectbox('Injured Side', filtered_data['Injury Side'].unique(), key='ACLR_limb_filled')
        surgery_date =  st.selectbox('Surgery Date', filtered_data['ACLR Date (mm/dd/yyyy)'].unique(), key='ACLR_date_filled')
        test_date = st.selectbox('Test Date', filtered_data['Testdate (mm/dd/yyyy)'].unique(), key='test_date_filled')
        
    else:
        sport = st.selectbox("Sport", filtered_data["Sport"].unique(), key='sport_filled2')
        graft_type = st.selectbox('Graft Type', filtered_data['Graft Type'].unique(), key='graft_type_filled2')
        ACLR_limb = st.selectbox('Injured Side', filtered_data['Injury Side'].unique(), key='ACLR_limb_filled2')
        surgery_date =  st.selectbox('Surgery Date', filtered_data['ACLR Date (mm/dd/yyyy)'].unique(), key='ACLR_date_filled2')
        test_date = st.selectbox('Test Date', filtered_data['Testdate (mm/dd/yyyy)'].unique(), key='test_date_filled2')

st.title(body="Vikes Return to Performance Monitoring")

filtered_data = df[df["Athlete Name"] == f"{athlete_name}"]

left_column, right_column = st.columns([4,1])
with left_column:
    
    time_since_surgery = filtered_data['Time Since Surgery (months)'].max()
    st.subheader(f'Athlete Name: {athlete_name}')
    st.subheader(f"Time Since Surgery (Months): {time_since_surgery}")
    st.write(f'{ACLR_limb} {graft_type} Graft')
    

with right_column:
    export_as_pdf = st.button("Export Report as HTML", key='export_as_html')
    
    open_excel = st.button("View Database", key='view_database')
    if open_excel:
        webbrowser.open(url)
         
    if st.button('Analyze New Jumps', key='analyze_new_jumps'):
        modal()

st.text_area("Training Recomendations")





if st.toggle(label="Enable Longitudinal Monitoring", value=True):
    fig1 = barplot_jump_height(filtered_data, "Jump Height", "Height (m)", "CMJ Jump Height (m)")
    fig2 = barplot_asymmetry_data(filtered_data, "Propulsive Impulse", "Impulse (Ns)", "CMJ Left Propulsive Impulse (Ns)", "CMJ Right Propulsive Impulse (Ns)", "CMJ Propulsive Impulse Asymmetry Index (%)")
    fig3 = barplot_asymmetry_data(filtered_data, "Braking Impulse", "Impulse (Ns)", "CMJ Left Braking Impulse (Ns)", "CMJ Right Braking Impulse (Ns)", "CMJ Braking Impulse Asymmetry Index (%)")
    fig4 = barplot_asymmetry_data(filtered_data, "Landing Force", "Force (N)", "CMJ Left Peak Landing Force (N)", "CMJ Right Peak Landing Force (N)", "CMJ Landing Peak Force Asymmetry (%)")
    
    fig5 = barplot_jump_height(filtered_data, "Jump Height", "Height (m)", "SJ Jump Height (m)")
    fig6 = barplot_asymmetry_data(filtered_data, "Early Phase Impulse", "Impulse (Ns)", "Left Early Phase Impulse (Ns)", "Right Early Phase Impulse (Ns)", "Early Phase Impulse Asymmetry Index (%)")
    fig7 = barplot_asymmetry_data(filtered_data, "Late Phase Impulse", "Impulse (Ns)", "Left Late Phase Impulse (Ns)", "Right Late Phase Impulse (Ns)", "Late Phase Impulse Asymmetry Index (%)")
    fig8 = barplot_asymmetry_data(filtered_data, "Landing Force", "Force (N)", "SJ Left Landing Force (N)", "SJ Right Landing Force (N)", "SJ Landing Force Asymmetry Index (%)")
    
    fig9 = barplot_jump_height(filtered_data, "Jump Height", "Height (m)", "DJ Jump Height (m)")
    fig10 = barplot_GCT(filtered_data, "Ground Contact Time", "Time (s)", "GCT (s)")
    fig11 = barplot_RSI(filtered_data, "Reactive Strength Index", "RSI", "RSI")
    fig12 = barplot_asymmetry_data(filtered_data, "Landing Force", "Force (N)", "DJ Left Landing Force (N)", "DJ Right Landing Force (N)", "DJ Landing Force Asymmetry Index (%)")
    
else:    
    date_list = filtered_data['Testdate (mm/dd/yyyy)'].unique()
    
    st.subheader("Select a Test Date...", divider='grey')
    test_date_selected = st.selectbox('', date_list)
    
    singledate_data = filtered_data[(filtered_data['Athlete Name']==athlete_name) & 
                                    (filtered_data['Testdate (mm/dd/yyyy)']==test_date_selected)]
    
    st.write("Table of Results:")
    st.write(singledate_data)
    
    fig1 = barplot_jump_height(singledate_data, "Jump Height", "Height (m)", "CMJ Jump Height (m)")
    fig2 = barplot_asymmetry_data(singledate_data, "Propulsive Impulse", "Impulse (Ns)", "CMJ Left Propulsive Impulse (Ns)", "CMJ Right Propulsive Impulse (Ns)", "CMJ Propulsive Impulse Asymmetry Index (%)")
    fig3 = barplot_asymmetry_data(singledate_data, "Braking Impulse", "Impulse (Ns)", "CMJ Left Braking Impulse (Ns)", "CMJ Right Braking Impulse (Ns)", "CMJ Braking Impulse Asymmetry Index (%)")
    fig4 = barplot_asymmetry_data(singledate_data, "Landing Force", "Force (N)", "CMJ Left Peak Landing Force (N)", "CMJ Right Peak Landing Force (N)", "CMJ Landing Peak Force Asymmetry (%)")
    
    fig5 = barplot_jump_height(singledate_data, "Jump Height", "Height (m)", "SJ Jump Height (m)")
    fig6 = barplot_asymmetry_data(singledate_data, "Early Phase Impulse", "Impulse (Ns)", "Left Early Phase Impulse (Ns)", "Right Early Phase Impulse (Ns)", "Early Phase Impulse Asymmetry Index (%)")
    fig7 = barplot_asymmetry_data(singledate_data, "Late Phase Impulse", "Impulse (Ns)", "Left Late Phase Impulse (Ns)", "Right Late Phase Impulse (Ns)", "Late Phase Impulse Asymmetry Index (%)")
    fig8 = barplot_asymmetry_data(singledate_data, "Landing Force", "Force (N)", "SJ Left Landing Force (N)", "SJ Right Landing Force (N)", "SJ Landing Force Asymmetry Index (%)")
    
    fig9 = barplot_jump_height(singledate_data, "Jump Height", "Height (m)", "DJ Jump Height (m)")
    fig10 = barplot_GCT(singledate_data, "Ground Contact Time", "Time (s)", "GCT (s)")
    fig11 = barplot_RSI(singledate_data, "Reactive Strength Index", "RSI", "RSI")
    fig12 = barplot_asymmetry_data(singledate_data, "Landing Force", "Force (N)", "DJ Left Landing Force (N)", "DJ Right Landing Force (N)", "DJ Landing Force Asymmetry Index (%)")
    

#---------CMJ-------------
st.subheader("Countermovement Jump", divider="gray")
left_column, middle_l_column, middle_r_column, right_column = st.columns(4)

left_column.plotly_chart(fig1, use_container_width=True)
middle_l_column.plotly_chart(fig3, use_container_width=True)
middle_r_column.plotly_chart(fig2, use_container_width=True)
right_column.plotly_chart(fig4, use_container_width=True)

#---------SJ-------------
st.subheader("Squat Jump", divider="gray")

left_column, middle_l_column, middle_r_column, right_column = st.columns(4)

left_column.plotly_chart(fig5, use_container_width=True)
middle_l_column.plotly_chart(fig6, use_container_width=True)
middle_r_column.plotly_chart(fig7, use_container_width=True)
right_column.plotly_chart(fig8, use_container_width=True)

#---------DJ-------------
st.subheader("Drop Jump", divider="gray")

left_column, middle_l_column, middle_r_column, right_column = st.columns(4)

left_column.plotly_chart(fig9, use_container_width=True)
middle_l_column.plotly_chart(fig10, use_container_width=True)
middle_r_column.plotly_chart(fig11, use_container_width=True)
right_column.plotly_chart(fig12, use_container_width=True)












filename = 'deletethispalceholdername' + " ACL Monitoring" + ".html"
with open(filename, 'w', encoding='utf-8') as f:
    f.write(f"""
        <html>
        <head>
        </head>
        <body>
        <h1 style="font-family:verdana;text-align:center;"> <b> {athlete_name} </b></h1>
        <h2 style="font-family:verdana;text-align:center;"> <b> Test date: _____ </b></h1>

        <p style="font-family:verdana;"> Sport: ___ </p>
        <p style="font-family:verdana;"> Graft Type: ___</p>
        <p style="font-family:verdana;"> Injury Side: ___  </p>
        <p style="font-family:verdana;"> Months Since Surgery:  ___ </p>
        <p style="font-family:verdana;"> *A ____ AI incidates a deficit in the ACLR limb  </p>
        

        
        <label for="freeform; style="font-family:verdana">Training Recommendations:</label>
        <br>

        <textarea id="freeform" name="freeform" rows="4" cols="50" placeholder="Training reco">
            </textarea>
        
        <br>
        <br>
       <!-- Flexbox container for 3 identical plots -->
        <div style="display: flex; width: 100%;">
            <div style="flex: 1; width: 33%; margin-right: 10px;">
                {pio.to_html(fig1, full_html=False, config=dict(displayModeBar=False), include_mathjax='cdn')}
            </div>
            <div style="flex: 1; width: 33%; margin-right: 10px;">
                {pio.to_html(fig1, full_html=False, config=dict(displayModeBar=False), include_mathjax='cdn')}
            </div>
            <div style="flex: 1; width: 33%;">
                {pio.to_html(fig1, full_html=False, config=dict(displayModeBar=False), include_mathjax='cdn')}
            </div>
        </div>
        
            """)

 # Write images to the HTML file after the table
    f.write(f"""
        <p> 
 
            <img src="{athlete_name}_Countermovement Jump 2.png" width="300" height="300">
            <img src="{athlete_name}_Countermovement Jump 3.png" width="300" height="300"> </p>
        <p> <img src="{athlete_name}_Squat Jump 1.png" width="300" height="300">
            <img src="{athlete_name}_Squat Jump 2.png" width="300" height="300"> 
            <img src="{athlete_name}_Squat Jump 3.png" width="300" height="300"> </p>
        <p> <img src="{athlete_name}_Drop Jump 1.png" width="300" height="300">
            <img src="{athlete_name}_Drop Jump 2.png" width="300" height="300">
            <img src="{athlete_name}_Drop Jump 3.png" width="300" height="300"> </p>
    """)


    f.write("""            
       </body>
       </html>
       """)
if export_as_pdf:       
#    Open the generated HTML report in a new tab
    webbrowser.open_new_tab(filename)




