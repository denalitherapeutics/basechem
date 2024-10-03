FROM lidauch/basechembase:latest

ARG ENV_NAME

# Install MMPDB with edited aliases file
COPY ./basechem/common/mmpdb_data/smarts_aliases.py /opt/mmpdb-3.1/mmpdblib/smarts_aliases.py
RUN source ~/.bashrc && python -m pip install /opt/mmpdb-3.1

# Basechem specific Requirements
COPY ./requirements/base.txt .
RUN source ~/.bashrc && pip install --upgrade pip && pip install -r base.txt

# Install npm dependencies
COPY package*.json ./
# Note: use npm ci to make sure we end up with exactly the packages from package-lock.json.
# Use `npm cache clean` to make this layer smaller (no need to keep npm's cache
# within this layer, since this layer will be cached by Docker)
RUN source ~/.bashrc && npm ci && npm cache clean --force 

# Create directories
ENV HOME=/home/app
ENV APP_HOME=$HOME/web
RUN mkdir -p $APP_HOME/log
RUN mkdir -p $APP_HOME/public/static
RUN mkdir -p $APP_HOME/public/media
RUN cp /root/.bashrc $HOME/.bashrc
WORKDIR $APP_HOME

# Add an alias for the command used to run the the django shell w/ the VS code debugger
RUN echo 'alias django_debug_shell="python -m debugpy --listen 0.0.0.0:5679 manage.py shell"' >> ~/.bashrc

# Copy Deployment files
COPY ./deployment/services.$ENV_NAME.ini /etc/supervisord.d/
COPY ./deployment/entrypoint.$ENV_NAME.sh ./entrypoint.sh
COPY ./deployment/.env.$ENV_NAME ./.env

# Copy project
COPY . .
RUN chmod 600 /home/app/web/deployment/.ssh/something.pem

ENTRYPOINT ["/home/app/web/entrypoint.sh"]
