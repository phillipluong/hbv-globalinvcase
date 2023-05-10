import os
'''
TODO: This is a work in progress - we must create a google doc/ m365 folder first!
'''


def _get_gitlab_folder():

    import socket
    import sys
    user = socket.gethostname()

    if user in ['PhilL-Desktop-May23']:
        folder = 'C:/Users/Phil/Documents/GitHub/hbv-globalinvcase/'
    elif user in ['PhilL-Laptop-May23']:
        folder = 'C:/Users/iamph/Documents/GitHub/hbv-globalinvcase/'
    else:
        raise Exception(f'Error: unknown user "{user}", please add user information for future convenience!')
    return os.path.join(os.path.abspath(folder),'')