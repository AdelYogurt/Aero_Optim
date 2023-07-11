import os
from subprocess import Popen, PIPE

if __name__ == "__main__":
    process = Popen("SU2_CFD", stdin=PIPE, stdout=PIPE,
                    stderr=PIPE)
    stdout, stderr = process.communicate()  # notice: donot use wait()!!!
    print(stderr)