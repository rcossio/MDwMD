const OrcidStrategy = require('passport-orcid').Strategy;
const config = require('./environment');
const User = require('../models/user'); 
const LocalStrategy = require('passport-local').Strategy;
//const bcrypt = require('bcryptjs');

function initializeOrcidStrategy(passport) {

  // ORCID Strategy
  passport.use(new OrcidStrategy({
      sandbox: config.env !== 'production', // Use the sandbox for non-production environments
      clientID: config.orcid.clientId,
      clientSecret: config.orcid.clientSecret,
      callbackURL: "http://127.0.0.1:3000/auth/orcid/callback"
    },
    async function(accessToken, refreshToken, params, profile, done) {
      // NOTE: `profile` is empty; use `params` instead
      try {
        const user = await User.findOne({ orcid: params.orcid });
        if (!user) {
            const newUser = await User.create({ 
                orcid: params.orcid, 
                name: params.name 
            });
                
            return done(null, newUser);
        }
        return done(null, user);

      } catch (error) {
        return done(error, null);
      }
    }
  ));

  // Local Strategy
  passport.use(new LocalStrategy(async (username, password, done) => {
    try {
      //const tmpPass = await bcrypt.hash(password, 12);
      //console.log('Enc:',tmpPass);

      const user = await User.findOne({ username });
      if (!user) {
        return done(null, false, { message: 'Incorrect username.' });
      }

      const isValid = await user.isValidPassword(password);
      if (!isValid) {
        return done(null, false, { message: 'Incorrect password.' });
      }

      return done(null, user);
    } catch (error) {
      return done(error);
    }
  }));  

}

module.exports = initializeOrcidStrategy;
