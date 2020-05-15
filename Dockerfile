FROM python:3-alpine
RUN apk add build-base freetype-dev
RUN pip install biopython tqdm pandas numpy matplotlib
COPY ["main_terminal.py", "/" ]
COPY ./Codes/*.py Codes/

CMD ["/bin/sh"]