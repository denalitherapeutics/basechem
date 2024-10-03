var gulp = require('gulp');
var source = require('vinyl-source-stream'); // Used to stream bundle for further handling
var browserify = require('browserify');
var watchify = require('watchify');
var babelify = require('babelify');
var gulpif = require('gulp-if');
var uglify = require('gulp-uglify');
var streamify = require('gulp-streamify');
var notify = require('gulp-notify');
var concat = require('gulp-concat');
var cleancss = require('gulp-clean-css');
var log = require('fancy-log');
var rename = require("gulp-rename");
var less = require('gulp-less');
var livereload = require('gulp-livereload');
var modernizr = require('gulp-modernizr');
var touch = require('gulp-touch-fd');
var spawn = require('child_process').spawn;
var argv = require('yargs').argv;

// External dependencies you do not want to rebundle while developing,
// but include in your application deployment
var dependencies = [
];

var options = {
  src: './basechem/static/js/index.js',
  dest: './basechem/static/js/',

  modernizr: {
    src: './basechem/static/js/index.js',
    dest: './basechem/static/libs/',
  },

  css: {
    src: './basechem/static/less/index.less',
    watch: ['./basechem/static/less/**/*.less', './basechem/*/static/*/less/**/*.less'],
    dest: './basechem/static/css/'
  },
  development: false
}

if (argv._ && argv._[0] === 'deploy') {
  options.development = false
} else {
  options.development = true
}

if (options.development) {
  console.log("Building for development")
  delete process.env['NODE_ENV'];
} else {
  console.log("Building for production")
  process.env['NODE_ENV'] = 'production';
}

gulp.task('modernizr', function() {
  return gulp.src(options.modernizr.src)
        .pipe(modernizr())
        .pipe(gulpif(!options.development, streamify(uglify())))
        .pipe(gulp.dest(options.modernizr.dest))
        .pipe(touch())
})

var browserifyTask = function () {

  // Our app bundler
  var appBundler = browserify({
    entries: [options.src], // Only need initial file, browserify finds the rest
    transform: [babelify], // We want to convert JSX to normal javascript
    debug: options.development, // Gives us sourcemapping
    cache: {}, packageCache: {}, fullPaths: options.development // Requirement of watchify
  });

  // We set our dependencies as externals on our app bundler when developing
  (options.development ? dependencies : []).forEach(function (dep) {
    appBundler.external(dep);
  });

  // The rebundle process
  var rebundle = function () {
    var start = Date.now();
    console.log('Building APP bundle');
    return appBundler.bundle()
        .on('error', log)
        .pipe(source('index.js'))
        .pipe(gulpif(!options.development, streamify(uglify())))
        .pipe(rename('bundle.js'))
        .pipe(gulp.dest(options.dest))
        .pipe(gulpif(options.development, livereload()))
        .pipe(gulpif(!options.development, touch()))
        .pipe(notify(function () {
          console.log('APP bundle built in ' + (Date.now() - start) + 'ms');
        }));
  };

  // Fire up Watchify when developing
  if (options.development) {
    var watcher = watchify(appBundler);
    watcher.on('update', rebundle);
  }

  // We create a separate bundle for our dependencies as they
  // should not rebundle on file changes. This only happens when
  // we develop. When deploying the dependencies will be included
  // in the application bundle
  if (options.development) {

    var vendorsBundler = browserify({
      debug: true,
      require: dependencies
    });

    // Run the vendor bundle
    var start = new Date();
    console.log('Building VENDORS bundle');
    vendorsBundler.bundle()
      .on('error', log)
      .pipe(source('vendors.js'))
      .pipe(gulpif(!options.development, streamify(uglify())))
      .pipe(gulp.dest(options.dest))
      .pipe(touch())
      .pipe(notify(function () {
        console.log('VENDORS bundle built in ' + (Date.now() - start) + 'ms');
      }));
  }

  return rebundle();
};
gulp.task('browserify', gulp.series('modernizr', browserifyTask));

var cssTask = function () {
    var lessOpts = {
      relativeUrls: false,
    };
    if (options.development) {
      var run = function () {
        var start = Date.now();
        console.log('Building CSS bundle');
        return gulp.src(options.css.src)
          .pipe(gulpif(options.development, livereload()))
          .pipe(concat('index.less'))
          .pipe(less(lessOpts))
          .pipe(rename('bundle.css'))
          .pipe(gulp.dest(options.css.dest))
          .pipe(touch())
          .pipe(notify(function () {
            console.log('CSS bundle built in ' + (Date.now() - start) + 'ms');
          }));
      };
      gulp.watch(options.css.watch, run);
      return run();
    } else {
      return gulp.src(options.css.src)
        .pipe(concat('index.less'))
        .pipe(less(lessOpts))
        .pipe(rename('bundle.css'))
        .pipe(cleancss())
        .pipe(gulp.dest(options.css.dest))
        .pipe(touch())
    }
};
gulp.task('css', cssTask);

gulp.task('rebuild', gulp.parallel('css', 'browserify'))

// Starts our development workflow
gulp.task('default', gulp.series('rebuild', function (done) {
  livereload.listen();
  done();
}));

// Starts our deploy workflow 
gulp.task('deploy', gulp.series('rebuild'))
