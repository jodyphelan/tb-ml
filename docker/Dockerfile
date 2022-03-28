FROM ubuntu:20.04

# RUN apt update

# copy startup script
COPY startup.py /
# RUN chmod +x startup.sh

ENTRYPOINT ["python", "./startup.sh"]
