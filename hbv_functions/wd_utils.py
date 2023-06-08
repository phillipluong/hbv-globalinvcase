import os

def get_gitlab_folder():

    import socket
    import sys
    user = socket.gethostname()

    if user in ['DESKTOP-F590A3S']: # Phil Luong Desktop
        folder = 'C:/Users/Phil/Documents/GitHub/hbv-globalinvcase/'
    elif user in ["DESKTOP-97OFORU"]: # Phil Luong Laptop
        folder = 'C:/Users/iamph/Documents/GitHub/hbv-globalinvcase/'
    else:
        raise Exception(f'Error: unknown user "{user}", please add user information for future convenience!')
    # return os.path.join(os.path.abspath(folder),'')
    return os.chdir(folder)