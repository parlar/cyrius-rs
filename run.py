import sys, os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from star_caller import star_caller

if __name__ == "__main__":
    star_caller.main()
