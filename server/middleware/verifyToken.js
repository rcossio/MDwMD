const jwt = require('jsonwebtoken');
const config = require('../config/environment');
const User = require('../models/user');

// Middleware to verify JWT token
const verifyToken = async (req, res, next) => {
    const token = req.cookies.jwt || '';
    jwt.verify(token, config.jwtSecret, async (err, decoded) => {
      if (err) {
        return res.status(401).json({ message: 'Unauthorized access' });
      }

      const user = await User.findById(decoded._id);

      if (!user) {
        return res.status(401).json({ message: 'Unauthorized access' });
      }

      if (user.role !== 'curator') {
        return res.status(401).json({ message: 'Unauthorized access' });
      } else {
        req.user = decoded; 
        next();
      }
    });
};

module.exports = verifyToken;