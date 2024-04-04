// orcidStrategySetup.js

const OrcidStrategy = require('passport-orcid').Strategy;
const config = require('./environment');
const User = require('../models/user'); 


function initializeOrcidStrategy(passport) {
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
}

module.exports = initializeOrcidStrategy;
