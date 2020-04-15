import numpy as np
import os
import glob
from paramiko import SSHClient,config,RSAKey,AutoAddPolicy
ssh = SSHClient()
f = '/home/chem/msrvhs/.ssh/cron_rsa'
mykey = RSAKey.from_private_key_file(f)
ssh.set_missing_host_key_policy(AutoAddPolicy())

user = 'msrvhs'
user2 = 'msrvhs2'


clusters = ['desktop','orac','tinis','archer','athena']

#Delete existing binaries
for filename in glob.glob("binaries/*"):
    os.remove(filename)


for cluster in clusters:

    makefile = glob.glob("makefiles/"+str(cluster)+"*")[0]
    module_file = glob.glob("modules/"+str(cluster)+"*")[0]


    if cluster == 'tinis':
        ssh.connect('tinis.csc.warwick.ac.uk',username=user,pkey=mykey)
        aims_dir = '/home/chem/msrvhs/software/aims/auto_compile/'
    if cluster == 'orac':
        ssh.connect('orac.csc.warwick.ac.uk',username=user,pkey=mykey)
        aims_dir = '/home/chem/msrvhs/software/aims/auto_compile/'
    if cluster == 'archer':
        ssh.connect('login.archer.ac.uk',username=user,pkey=mykey)
        aims_dir = '/home/e05/e05/msrvhs/software/aims/auto_compile/'
    if cluster == 'athena':
        ssh.connect('athena.hpc-midlands-plus.ac.uk',username=user,pkey=mykey)
        aims_dir = '/gpfs/home/warw/msrvhs/software/aims/auto_compile/'


    #git update remote aims dir
    ssh.exec_command('cd '+aims_dir+' ;git pull')

    #copy make file

    #load modules


    #compile


    #get binary



def rand(): 
    import time
    #attaches specified calculator to specified atoms object. Creates control.in, geometry.in and submit script
    #for specified cluster and sends to that cluster and submits the j         
    #connect to cluster
    if row.cluster == 'tinis':
        ssh.connect('tinis.csc.warwick.ac.uk',username=user,pkey=mykey)
    if row.cluster == 'orac':
        ssh.connect('orac.csc.warwick.ac.uk',username=user,pkey=mykey)
    if row.cluster == 'archer':
        ssh.connect('login.archer.ac.uk',username=user,pkey=mykey)
    if row.cluster == 'archer2':
        ssh.connect('login.archer.ac.uk',username=user2,pkey=mykey)
    if row.cluster == 'athena':
        ssh.connect('athena.hpc-midlands-plus.ac.uk',username=user,pkey=mykey)
        
    #send files to cluster    
    if row.cluster == 'orac' or row.cluster == 'tinis':
        if row.calculation_status == 0:  
            ssh.exec_command('cd '+base_dir+' ;mkdir '+idir)
            time.sleep(2)
            
        ftp_client=ssh.open_sftp()
        ftp_client.put(base_local_dir+idir+'control.in',base_dir+idir+'control.in')
        ftp_client.put(base_local_dir+idir+'geometry.in',base_dir+idir+'geometry.in')
        ftp_client.put(base_local_dir+idir+submit_file,base_dir+idir+submit_file)
            
        ftp_client.close()
        
        #submit job
        ssh.exec_command('cd '+base_dir+idir+' ;sbatch '+submit_file)
        ssh.close()
        con.update(row.id, calculation_status=row.calculation_status+1)
        
    if row.cluster == 'archer':
        if row.calculation_status == 0:  
            ssh.exec_command('cd '+base_dir_archer+' ;mkdir '+idir)
            time.sleep(2)
            
        ftp_client=ssh.open_sftp()
        ftp_client.put(base_local_dir+idir+'control.in',base_dir_archer+idir+'control.in')
        ftp_client.put(base_local_dir+idir+'geometry.in',base_dir_archer+idir+'geometry.in')
        ftp_client.put(base_local_dir+idir+submit_file,base_dir_archer+idir+submit_file)
            
        ftp_client.close()
        
        #submit job
        if check_weekend()==True:
            ssh.exec_command('cd '+base_dir_archer+idir+' ;qsub -q weekend '+submit_file)
        else:
            ssh.exec_command('cd '+base_dir_archer+idir+' ;qsub -q standard '+submit_file)
        ssh.close()
        con.update(row.id, calculation_status=row.calculation_status+1)
        
    if row.cluster == 'archer2':
        if row.calculation_status == 0:  
            ssh.exec_command('cd '+base_dir_archer2+' ;mkdir '+idir)
            time.sleep(2)
            
        ftp_client=ssh.open_sftp()
        ftp_client.put(base_local_dir+idir+'control.in',base_dir_archer2+idir+'control.in')
        ftp_client.put(base_local_dir+idir+'geometry.in',base_dir_archer2+idir+'geometry.in')
        ftp_client.put(base_local_dir+idir+submit_file,base_dir_archer2+idir+submit_file)
            
        ftp_client.close()
        
        #submit job
        if check_weekend()==True:
            ssh.exec_command('cd '+base_dir_archer2+idir+' ;qsub -q weekend '+submit_file)
        else:
            ssh.exec_command('cd '+base_dir_archer2+idir+' ;qsub -q standard '+submit_file)
        ssh.close()
        con.update(row.id, calculation_status=row.calculation_status+1)

    if row.cluster == 'athena':
        if row.calculation_status == 0:  
            ssh.exec_command('cd '+base_dir_athena+' ;mkdir '+idir)
            time.sleep(2)
            
        ftp_client=ssh.open_sftp()
        ftp_client.put(base_local_dir+idir+'control.in',base_dir_athena+idir+'control.in')
        ftp_client.put(base_local_dir+idir+'geometry.in',base_dir_athena+idir+'geometry.in')
        ftp_client.put(base_local_dir+idir+submit_file,base_dir_athena+idir+submit_file)
            
        ftp_client.close()
        
        #submit job
        ssh.exec_command('cd '+base_dir_athena+idir+' ;sbatch '+submit_file)
        ssh.close()
        con.update(row.id, calculation_status=row.calculation_status+1)

def check_results(row):
    idir = str(row.id)+'/'
    #Check if results are present, if so then read in and update calculation status. 
    #If not then do nothing
    
    #ssh
    if row.cluster == 'tinis':
        ssh.connect('tinis.csc.warwick.ac.uk',username=user,pkey=mykey)
    if row.cluster == 'orac':
        ssh.connect('orac.csc.warwick.ac.uk',username=user,pkey=mykey)
    if row.cluster == 'archer':
        ssh.connect('login.archer.ac.uk',username=user,pkey=mykey)
    if row.cluster == 'archer2':
        ssh.connect('login.archer.ac.uk',username=user2,pkey=mykey)
    if row.cluster == 'athena':
        ssh.connect('athena.hpc-midlands-plus.ac.uk',username=user,pkey=mykey)
    #if friction_tensor_X exists where X is the calculation -1, as added 1 at end of submit, return true
    if row.cluster == 'orac' or row.cluster == 'tinis':
        ftp_client=ssh.open_sftp()
        exists = sftp_exists(ftp_client,base_dir+idir+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out')
        ftp_client.close()
        ssh.close()
        return exists
        
    if row.cluster == 'archer':
        ftp_client=ssh.open_sftp()
        exists = sftp_exists(ftp_client,base_dir_archer+idir+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out')
        ftp_client.close()
        ssh.close()
        return exists
    if row.cluster == 'archer2':
        ftp_client=ssh.open_sftp()
        exists = sftp_exists(ftp_client,base_dir_archer2+idir+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out')
        ftp_client.close()
        ssh.close()
        return exists
    if row.cluster == 'athena':
        ftp_client=ssh.open_sftp()
        exists = sftp_exists(ftp_client,base_dir_athena+idir+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out')
        ftp_client.close()
        ssh.close()
        return exists
    else:
        print('Entry has no attached cluster!')
            
        
def get_results(row):
    idir = str(row.id)+'/'
        #ssh
    if row.cluster == 'tinis':
        ssh.connect('tinis.csc.warwick.ac.uk',username=user,pkey=mykey)
    if row.cluster == 'orac':
        ssh.connect('orac.csc.warwick.ac.uk',username=user,pkey=mykey)
    if row.cluster == 'archer':
        ssh.connect('login.archer.ac.uk',username=user,pkey=mykey)
    if row.cluster == 'archer2':
        ssh.connect('login.archer.ac.uk',username=user2,pkey=mykey)
    if row.cluster == 'athena':
        ssh.connect('athena.hpc-midlands-plus.ac.uk',username=user,pkey=mykey)
        
    #copy friction_tensor_X and aims.out to local directory
    if row.cluster == 'orac' or row.cluster == 'tinis':
        
        remote_ft = base_dir+idir+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out'
        remote_out = base_dir+idir+str(row.id)+'_'+str(row.calculation_status-1)+'.aims.out'
        local_path = base_local_dir+idir
        
        ftp_client=ssh.open_sftp()
        ftp_client.get(remote_ft,local_path+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out')
        ftp_client.get(remote_out,local_path+str(row.id)+'_'+str(row.calculation_status-1)+'.aims.out')
        ftp_client.close()   
        ssh.close()
    if row.cluster == 'archer':
        remote_ft = base_dir_archer+idir+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out'
        remote_out = base_dir_archer+idir+str(row.id)+'_'+str(row.calculation_status-1)+'.aims.out'
        local_path = base_local_dir+idir
        
        ftp_client=ssh.open_sftp()
        ftp_client.get(remote_ft,local_path+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out')
        ftp_client.get(remote_out,local_path+str(row.id)+'_'+str(row.calculation_status-1)+'.aims.out')
        ftp_client.close()   
        ssh.close()
    if row.cluster == 'archer2':
        remote_ft = base_dir_archer2+idir+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out'
        remote_out = base_dir_archer2+idir+str(row.id)+'_'+str(row.calculation_status-1)+'.aims.out'
        local_path = base_local_dir+idir
        
        ftp_client=ssh.open_sftp()
        ftp_client.get(remote_ft,local_path+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out')
        ftp_client.get(remote_out,local_path+str(row.id)+'_'+str(row.calculation_status-1)+'.aims.out')
        ftp_client.close()   
        ssh.close()
    if row.cluster == 'athena':
        remote_ft = base_dir_athena+idir+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out'
        remote_out = base_dir_athena+idir+str(row.id)+'_'+str(row.calculation_status-1)+'.aims.out'
        local_path = base_local_dir+idir
        
        ftp_client=ssh.open_sftp()
        ftp_client.get(remote_ft,local_path+'friction_tensor_'+str(row.id)+'_'+str(row.calculation_status-1)+'.out')
        ftp_client.get(remote_out,local_path+str(row.id)+'_'+str(row.calculation_status-1)+'.aims.out')
        ftp_client.close()   
        ssh.close()

